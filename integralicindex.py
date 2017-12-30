from __future__ import print_function

import tempfile
import os
import argparse
import time

import astropy.io.fits as pyfits
import pilton
from ddosa import remove_withtemplate
import numpy as np
from collections import defaultdict
import timesystem


class ICTree(object):
    def __init__(self,suffix=""):
        self.icroot=os.environ['CURRENT_IC']
        self.suffix=suffix

    #@property
    def get_ibisicroot(self,DS):
        if DS=="ISGR-EBDS-MOD":
            return self.icroot+"/ic/ibis/rsp"
        
        if DS.endswith("RSP"):
            return self.icroot+"/ic/ibis/rsp"

        if DS.endswith("MOD"):
            return self.icroot+"/ic/ibis/mod"

    @property
    def idxicroot(self):
        return self.icroot+"/idx/ic/"

    @property
    def icmaster(self):
        return self.get_icmaster(self.suffix)
    
    def get_icmaster(self,version=None):
        return self.icroot+"/idx/ic/ic_master_file"+("_"+version if version is not None and version != "" else "")+".fits"

    def DS_to_fn_prefix(self,DS):
        return DS.lower().replace("-","_").replace(".","")
    
    def DS_to_mnemcol(self,DS):
        return DS.replace("-","_").replace(".","")

    def DS_to_fn(self,DS,serial=0):
        return self.get_ibisicroot(DS)+"/"+self.DS_to_fn_prefix(DS)+"_%.4i.fits"%serial

    def DS_to_idx_fn(self,DS):
        if self.suffix=="":
            return self.idxicroot+"/"+DS+"-IDX.fits"
        else:
            return self.idxicroot+"/"+DS+"-IDX_%s.fits"%self.suffix


    def create_index_empty(self,DS):
        dc=pilton.heatool("dal_create")
        dc["obj_name"]=self.DS_to_idx_fn(DS)
        dc["template"]=DS+"-IDX.tpl"
        remove_withtemplate(dc["obj_name"].value+"("+dc["template"].value+")")
        dc.run()
        
    def get_file_DS(self,fn):
        f=pyfits.open(fn)
        if len(f)>2:
            print("we refuse to deal with indexed IC ds, as they increase the amount of suffering in the world")
            raise Exception("too many extensions %i %s"%(len(f),repr(f)))
        if len(f)<2:
            print("")
            raise Exception("too few extensions %i %s"%(len(f),repr(f)))
        return f[1].header['EXTNAME']


    def attach_ds(self,fn,serial=0):
        DS=self.get_file_DS(fn)
        f_ds=pyfits.open(fn)
        f_ds.writeto(self.DS_to_fn(DS,serial),clobber=True)

        da=pilton.heatool("dal_attach")
        da['Parent']=self.DS_to_idx_fn(DS)
        da['Child1']=self.DS_to_fn(DS,serial)
        da.run()

        f_idx=pyfits.open(da['Parent'].value)

        f_ds[1].header['VSTOP']=99999

        for k in ['VERSION','VSTART','VSTOP']:
            f_idx[1].data[-1][k]=f_ds[1].header[k]
            print(k,f_idx[1].data[-1][k],f_ds[1].header[k])

        f_idx.writeto(da['Parent'].value,clobber=True)

        dv=pilton.heatool("dal_verify")
        dv["indol"]=self.DS_to_idx_fn(DS)
        dv['checksums']="yes"
        dv['backpointers']="yes"
        dv['detachother']="yes"
        dv.run()



    def init_icmaster(self):
        f=pyfits.open(self.get_icmaster("osa102"))
        
        print("columns was:",len(f[3].columns))

        #print([c.name for c in f[3].columns])

        nhdu=pyfits.BinTableHDU.from_columns(f[3].columns+pyfits.ColDefs([
            pyfits.Column(self.DS_to_mnemcol(DS),"1I") for DS in self.icstructures.keys()
             if self.DS_to_mnemcol(DS) not in [c.name for c in f[3].columns]]))

        for k,v in f[3].header.items():
            if k.startswith("TFORM"): continue
            if k.startswith("TTYPE"): continue
            if k.startswith("TFIELDS"): continue
            if k.startswith("NAXIS"): continue
            #print ":",k,"x",v
            nhdu.header[k]=v

        for DS,icstructure in self.icstructures.items():
            versions=sorted(set([icfile['version'] for icfile in icstructure]))
            if len(versions)!=1:
                raise Exception("inconsistent versions for "+DS+" "+repr(versions))
            nhdu.data[0][self.DS_to_mnemcol(DS)]=versions[0]

        f[3]=nhdu

        print("columns now:",len(f[3].columns))

        f.writeto(self.get_icmaster(self.suffix),clobber=True)

    def create_index_from_list(self,DS,fns=None,fns_list=None,update=True,recreate=True):
        if fns_list is None:
            if fns is None:
                raise Exception("what?")

            fns_list_handle,fns_list_fn=tempfile.mkstemp()
            print("\n".join(fns))
            os.write(fns_list_handle,"\n".join(fns)+"\n")

        else:
            if fns is not None:
                raise Exception("what?")

        print("files:",fns)

        if os.path.exists(self.DS_to_idx_fn(DS)) and recreate:
            os.remove(self.DS_to_idx_fn(DS))

        da=pilton.heatool("txt2idx")
        da['index']=self.DS_to_idx_fn(DS)
        da['template']=DS+"-IDX.tpl"
        da['update']=1 if update else 0
        da['element']=fns_list_fn
        da.run()

        f=pyfits.open(da['index'].value)
        f[1].header['CREATOR']="Volodymyr Savchenko"
        f[1].header['CONFIGUR']="dev"
        f.writeto(da['index'].value,clobber=True)

    def attach_idx_to_master(self,DS):
        da=pilton.heatool("dal_attach")
        da['Parent']=self.icmaster+"[2]"
        da['Child1']=self.DS_to_idx_fn(DS)
        da.run()

    icstructures=defaultdict(list)

    def get_icfile_validity_rev(self,f,unique=True,first=True,middle=False):
        vstart=self.find_key(f,"VSTART")
        vstop=self.find_key(f,"VSTOP")

        rev_start=int(timesystem.converttime("IJD",vstart,"REVNUM"))

        if first:
            return rev_start
        
        rev_stop=int(timesystem.converttime("IJD",vstop,"REVNUM"))

        if middle:
            return round(0.5*(rev_start+rev_stop))
    
        if not unique:
            return (rev_start,rev_stop)

        if rev_start != rev_stop:
            raise Exception("this file has many revs",rev_start,rev_stop)

        return rev_start
        
    def find_version(self,f):
        return self.find_key(f,"VERSION")
        
    def find_key(self,f,k,unique=True):
        values=[]
        for e in f:
            if k in e.header:
                values.append(e.header[k])

        if not unique:
            return values

        values_unique=sorted(set(values))
        if len(values_unique)>1:
            raise Exception("this file has many versions",values_unique)
        if len(values_unique)==0:
            raise Exception("this file has no versions")
        return values_unique[0]

    def add_icfile(self,icfile):
        print("requested to add",icfile)
        DS=self.get_file_DS(icfile)
        print(icfile,"as",DS)

        hash_fn=os.path.dirname(os.path.abspath(icfile))+"/hash.txt"
        if os.path.exists(hash_fn):
            hashe=open(hash_fn).read()
        else:
            hashe=""

        f=pyfits.open(icfile)
        rev=self.get_icfile_validity_rev(f)

        if rev<0 or rev>9000: # over 9000!!
            serial=1
        else:
            serial=rev

        self.icstructures[DS].append(dict(
                    origin_filename=icfile,
                    version=self.find_version(f),
                    serial=serial,
                    hashe=hashe,
                    ))

    def write(self):
        self.init_icmaster()

        for DS,icfiles in self.icstructures.items():
            print(DS)

            filelist=[]
            for icfile in icfiles:
                print("IC file",icfile)
                ic_store_filename=self.DS_to_fn(DS,serial=icfile['serial'])
                print("store in IC as ",ic_store_filename)

                f_ds=pyfits.open(icfile['origin_filename'])
                f_ds[1].header['VSTOP']=99999
                f_ds.writeto(ic_store_filename,clobber=True)

                version_store=os.path.dirname(os.path.abspath(ic_store_filename))+"/.version."+os.path.basename(ic_store_filename)
                print("version store",version_store)
                open(version_store,"w").write(icfile['hashe'])
                icfile['size']=os.path.getsize(ic_store_filename)
                icfile['ic_store_filename']=ic_store_filename
                icfile['version_store']=version_store
    
                filelist.append(ic_store_filename)

            print("file list",filelist)

            idx_fn=self.DS_to_idx_fn(DS)
            if os.path.exists(idx_fn):
                print("index exists:",idx_fn,"OVERWRITING")

            self.create_index_from_list(DS,fns=filelist)
            self.attach_idx_to_master(DS)

            self.write_version()

    def write_version(self):
        open(self.icroot+"/idx/ic/version","w").write(time.strftime("%Y-%m-%dT%H:%M:%S"))
        open(self.icroot+"/ic/ibis/version","w").write(time.strftime("%Y-%m-%dT%H:%M:%S"))

    def summarize(self):
        for DS,icfiles in self.icstructures.items():
            print(DS,len(icfiles),"%.5lg"%(sum([k['size'] for k in icfiles])/1024./1024.),"Mb")



def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('icfile', metavar='new IC file',nargs="*",
            help='new IC file')
    parser.add_argument('-f','--file', metavar='new IC file', action='append',
            help='new IC file')
    parser.add_argument('-s','--suffix', metavar='suffix',
            help='new IC file',default="")
    #parser.add_argument('clobber-index', action='store_true',
    #        help='clobber index')

    args = parser.parse_args()

    ictree=ICTree(args.suffix)
    ictree.set_suffix=args.suffix

    for fn in args.file:
        print("from",fn)
        with open(fn) as f:
            for icfile in f:
                try:
                    ictree.add_icfile(icfile.strip())
                except Exception as e:
                    print("failed to add",icfile.strip(),e)

    for icfile in args.icfile:
        try:
            ictree.add_icfile(icfile)
        except Exception as e:
            print("failed to add",icfile)

    ictree.write()
    ictree.summarize()

 #       ictree.create_index_empty(DS)
 #       ictree.attach_ds(icfile)
 #       ictree.attach_to_master(DS)

if __name__ == "__main__":
    main()
