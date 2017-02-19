import httplib2
import urllib, logging
import time
import json


def authenticate():
    # Get an instance of a logger
    #logger = logging.getLogger(__name__)
    return newt_request('auth', 'POST', {'username' : 'pdacs', 'password' : 'e+oeJH4Hrtj4CMYee'}, )

def execute_job(pbsPath, cookieStr):
    return newt_request('queue/carver', 'POST', {'jobfile' : pbsPath}, cookieStr)


def newt_request(url, req_method, params=None, cookie_str=None):
    
    newt_base_url='https://newt.nersc.gov/newt/'
    
    full_url = newt_base_url+url
    conn = httplib2.Http()#disable_ssl_certificate_validation=True)
    
    # Massage inputs
    if cookie_str:
        headers={'Cookie': cookie_str}
    else:
        headers=None
    
    if type(params) is dict:
        body=urllib.urlencode(params)
    elif (type(params) is str) or (type(params) is unicode):
        body=params
    else:
        body=None
    #logger.debug("NEWT: %s %s"%(req_method,full_url))
    
    response, content = conn.request(full_url, req_method, body=body, headers=headers)

    #logger.debug("NEWT response: %s"%response.status)
    
    return (response, content) 
    
def getJobID(conn):
    response = json.loads("[" + conn + "]")
    return response[0]['jobid']
    
def getJobStatus(jobid, cookieStr):
    resStat, conStat = newt_request('queue/carver/' + jobid, 'GET', None, cookieStr)
    response = json.loads("[" + conStat + "]")
    response[0]['status']
    return response[0]['status']
    
def waitToFinish(jobid, cookieStr):
    status = getJobStatus(jobid, cookieStr)
    while status in "NBRQH":
        time.sleep(5)
        status = getJobStatus(jobid, cookieStr)
    return status
