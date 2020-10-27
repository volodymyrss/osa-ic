import subprocess
import os
import sys

target=sys.argv[1]

isdc_host="isdc-nx00.isdc.unige.ch"
isdc_root="/home/isdc/savchenk/osa11_deployment/ddcache"
local_root=os.environ['INTEGRAL_DDCACHE_ROOT']


subprocess.check_call([
                    "ssh",
                    isdc_host,
                    "mkdir -p "+isdc_root+"/"+target
                    ])

subprocess.check_call([
                    "rsync",
                    "-avu",
                    local_root+"/"+target+"/",
                    isdc_host+":"+isdc_root+"/"+target,
                    ])
