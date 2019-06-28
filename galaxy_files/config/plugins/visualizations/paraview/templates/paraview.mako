<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>${hda.name} | ${visualization_name}</title>

## ----------------------------------------------------------------------------
<link type="text/css" rel="stylesheet" media="screen" href="${h.url_for('/galaxy-pdacs-dev/static/scripts/libs/lib_nersc_pvweb/lib_nersc_pvweb_ui.css')}">
<script type="text/javascript" src="${h.url_for('/galaxy-pdacs-dev/static/scripts/libs/lib_nersc_pvweb/lib_nersc_pvweb.js')}"></script>
<script type="text/javascript" src="${h.url_for('/galaxy-pdacs-dev/static/scripts/libs/lib_nersc_pvweb/lib_nersc_pvweb_ui.js')}"></script>

## ----------------------------------------------------------------------------
<script type="text/javascript" src="https://newt.nersc.gov/js/jquery-1.7.2.min.js"></script>
<script type="text/javascript" src="https://newt.nersc.gov/js/newt.js"></script>
<script language="Javascript">
    $(document).ready( function() {
        // disable logging
        var console = {};
        console.log = function(){};
        window.console = console;

        // configure core
        lib_nersc_pvweb.set_web_host('portal-auth.nersc.gov')
        lib_nersc_pvweb.set_paraview_version('4.2.0-PDACS')
        lib_nersc_pvweb.set_paraview_prefix('/usr/common/graphics/ParaView')
        // configure ui
        lib_nersc_pvweb_ui.initialize()
        lib_nersc_pvweb_ui.set_web_app('https://portal-auth.nersc.gov:8000/apps/pvweb_visualizer.html')
        lib_nersc_pvweb_ui.create_submit_form($('#form_area'), '${hda.name}', '${hda.file_name}')
        lib_nersc_pvweb_ui.create_job_table($('#job_table'))
    });
</script>

</head>


## ----------------------------------------------------------------------------
<body>
<table border="0">
<tr>
<td id="form_area" style="vertical-align:top"></td>
<td id="job_table" style="vertical-align:top"></td>
</tr>
</table>
</script>    
</body>
