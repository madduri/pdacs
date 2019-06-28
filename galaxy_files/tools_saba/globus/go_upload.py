#!/usr/bin/env python

# Copyright 2010 University of Chicago
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Globus Online Upload Tool

Constructs a new Galaxy Dataset (Datasource?) from a GridFTP transfer.

python go_upload.py USERNAME -k ~/.globus/userkey.pem \
           -c ~/.globus/usercert.pem \
           -C ~/.globus/certificates/1c3f2ca8.0
           --source-ep=my#ep
           --source-path=/~/somefile.txt
           --deadline=5
"""

from datetime import datetime, timedelta
import sys
import time

from globusonline.transfer.api_client import SimpleTransfer
from globusonline.transfer.api_client import process_args, TransferAPIClient

# TransferAPIClient instance.
api = None
output = sys.stdout

def transfer(source_ep, source_filepath, destination_ep,
             destination_filepath, deadline):
    """d
    Do a bunch of API calls and display the results. Does a small transfer
    between tutorial endpoints, but otherwise does not modify user data.

    Uses module global API client instance.
    """
    display_activation(source_ep)
    display_activation(destination_ep)

    code, message, data = api.transfer_generate_id()
    transfer_id = data["value"]
    deadline = datetime.utcnow() + timedelta(minutes=10)
    t = SimpleTransfer(transfer_id, source_ep, destination_ep, deadline)
    t.add_item(source_filepath, destination_filepath)
    code, reason, data = api.transfer(t)
    task_id = data["task_id"]

    display_tasksummary(); print
    display_task(task_id); print

    # wait for the task to complete, and see the summary and lists
    # update
    if wait_for_task(task_id):
        print >>output,"=== After completion ==="
        display_tasksummary(); print
        display_task(task_id); print
        display_ls(destination_ep); print


def display_activation(endpoint_name):
    print >>output,"=== Endpoint pre-activation ==="
    display_endpoint(endpoint_name)
    print
    code, reason, result = api.endpoint_activate(endpoint_name, None,
                                                 if_expires_in=600)
    if result.code.startswith("AutoActivationFailed"):
        print >>output,"Auto activation failed, ls and transfers will likely fail!"
    print >>output,"result: %s (%s)" % (result.code, result.message)
    print >>output,"=== Endpoint post-activation ==="
    display_endpoint(endpoint_name)
    print


def display_tasksummary():
    code, reason, data = api.tasksummary()
    print >>output,"Task Summary for %s:" % api.username
    for k, v in data.iteritems():
        if k == "DATA_TYPE":
            continue
        print >>output,"%3d %s" % (int(v), k.upper().ljust(9))


def display_tasks(max_age=None):
    """
    @param max_age: only show tasks requested at or after now - max_age.
    @type max_age: timedelta
    """
    kwargs = {}
    if max_age:
        min_request_time = datetime.utcnow() - max_age
        # filter on request_time starting at min_request_time, with no
        # upper limit on request_time.
        kwargs["request_time"] = "%s," % min_request_time

    code, reason, tasks = api.tasks(**kwargs)
    print >>output,"Tasks for %s:" % api.username
    for task in tasks["DATA"]:
        print >>output,"Task %s:" % task["task_id"]
        _print_task(task)

def _print_task(data, indent_level=0):
    """
    Works for tasks and subtasks, since both have a task_id key
    and other key/values are printed by iterating through the items.
    """
    indent = " " * indent_level
    indent += " " * 2
    for k, v in data.iteritems():
        if k in ("DATA_TYPE", "LINKS"):
            continue
        print >>output,indent + "%s: %s" % (k, v)

def display_task(task_id, show_subtasks=True):
    code, reason, data = api.task(task_id)
    print >>output,"Task %s:" % task_id
    _print_task(data, 0)

    if show_subtasks:
        code, reason, data = api.subtasks(task_id)
        subtasks = data["DATA"]
        for t in subtasks:
            print >>output,"  subtask %s:" % t["task_id"]
            _print_task(t, 4)

def wait_for_task(task_id, timeout=120):
    status = "ACTIVE"
    while timeout and status == "ACTIVE":
        code, reason, data = api.task(task_id, fields="status")
        status = data["status"]
        time.sleep(1)
        timeout -= 1

    if status != "ACTIVE":
        print >>output,"Task %s complete!" % task_id
        return True
    else:
        print >>output,"Task still not complete after %d seconds" % timeout
        return False

def display_endpoint(endpoint_name):
    code, reason, data = api.endpoint(endpoint_name)
    _print_endpoint(data)

def display_ls(endpoint_name, path=""):
    code, reason, data = api.endpoint_ls(endpoint_name, path)
    # Server returns canonical path; "" maps to the users default path,
    # which is typically their home directory "/~/".
    path = data["path"]
    print >>output,"Contents of %s on %s:" % (path, endpoint_name)
    headers = "name, type, permissions, size, user, group, last_modified"
    headers_list = headers.split(", ")
    print >>output,headers
    for f in data["DATA"]:
        print >>output,", ".join([str(f[k]) for k in headers_list])


def _print_endpoint(ep):
    name = ep["canonical_name"]
    print >>output,name
    if ep["activated"]:
        print >>output,"  activated (expires: %s)" % ep["expire_time"]
    else:
        print >>output,"  not activated"
    if ep["public"]:
        print >>output,"  public"
    else:
        print >>output,"  not public"
    if ep["myproxy_server"]:
        print >>output,"  default myproxy server: %s" % ep["myproxy_server"]
    else:
        print >>output,"  no default myproxy server"
    servers = ep["DATA"]
    print >>output,"  servers:"
    for s in servers:
        print >>output,"    " + s["uri"],
        if s["subject"]:
            print >>output," (%s)" % s["subject"]
        else:
            print >>output, ""


def display_endpoints():
    code, reason, endpoints = api.endpoints(limit=100)
    print >>output, ("Found %d endpoints for user %s:" %
                     (endpoints["length"], api.username))
    for ep in endpoints["DATA"]:
        _print_endpoint(ep)


if __name__ == '__main__':
    options, args = process_args()
    api = TransferAPIClient(args[0], server_ca_file=options.server_ca_file,
                            cert_file=options.cert_file,
                            key_file=options.key_file,
                            saml_cookie=options.saml_cookie,
                            base_url=options.base_url)
    output = options.output
    display_endpoints()
    transfer(options.source_ep, options.source_path,
             options.destination_ep, options.destination_path,
             options.deadline)
