from __future__ import print_function

import dataanalysis as da
import ddosa

import os
from os import access, R_OK
from os.path import isfile
import re
import glob
import random            

import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy import units as u

import pilton

def converttime(informat,intime,outformat):
    ct=pilton.heatool("converttime")
    ct['informat']=informat
    ct['intime']=intime
    ct['outformat']=outformat
    ct.run()

    c=dict(re.findall("Output Time\((.*)\): (.*)",ct.output))
    if outformat=="":
        return c
    else:
        return c[outformat]

class TimeDirectionScWList(ddosa.DataAnalysis):
    coordinates=dict(RA=None,DEC=None,radius=None)
    timespan=dict(T1=None,T2=None)
    max_pointings=10
    randomize_pick=True

    allow_alias=True

    version="v1.1"

    scwversion="latest"

    def get_version(self):
        v=self.get_signature()+"."+self.version

        try:
            v+="c_%(RA).5lg_%(DEC).5lg_%(radius).5lg.t_%(T1)s_%(T2)s_maxpoint%(max_pointings)i"%(
                            dict(self.coordinates.items()+ \
                                 self.timespan.items()+ \
                                 dict(max_pointings=(self.max_pointings if self.max_pointings is not None else 1000)).items())
                        )
            v+="..."+self.timespan['T1'][:4]
            if self.timespan['T1'][:4]!=self.timespan['T2'][:4]:
                v+="_"+self.timespan['T2'][:4]
            else:
                v+="..."+self.timespan['T1'][5:7]
                if self.timespan['T1'][5:7]!=self.timespan['T2'][5:7]:
                    v+="_"+self.timespan['T2'][5:7]

        except Exception:
            print("failed coordinates",self.coordinates)
            print("failed timespan", self.timespan)
            print("failed max_pointings", self.max_pointings)
            v+=".UNSET" # TODO make it generic


        if self.randomize_pick:
            v+="...randompick"

        if self.scwversion!="any":
            v+="...scwversion."+self.scwversion

        return v


    def latest_index(self, rbp):
        index_fn = sorted(glob.glob(rbp+"/idx/scw/GNRL-SCWG-GRP-IDX_*"))[-1]

        index_version = re.search("/GNRL-SCWG-GRP-IDX_(.*?).fits(.gz)?", index_fn).groups()[1]

        print("searching for latest index in", rbp)
        print("found latest index", index_fn, index_version)

        return index_fn, index_version

    def scw_data(self, rbp, scwversion):
        index_fn, index_version = self.latest_index(rbp)


        scw_index = fits.open(index_fn)[1].data

        return self.extract_from_index(scw_index)

    def scw_data_cons(self):
        return self.scw_data(os.environ['INTEGRAL_DATA'], "001")

    def scw_data_nrt(self):
        return self.scw_data(os.environ['REP_BASE_PROD_NRT'], "000")
    
    def scw_data_latest(self):
        nrt_index_fn, nrt_index_version = self.latest_index(os.environ['REP_BASE_PROD_NRT'])
        cons_index_fn, cons_index_version = self.latest_index(os.environ['INTEGRAL_DATA'])

        nrt_scw_index = fits.open(nrt_index_fn)[1].data
        cons_scw_index = fits.open(cons_index_fn)[1].data

        print()

        raise Exception()

        #return self.scw_data(, "latest")

    def main(self):
        scw_cons=self.scw_data_cons()
        if self.scwversion=="001":
            print("instructed to use CONS")
            self.scwlistdata=scw_cons
            return
        
        if self.scwversion=="000":
            print("instructed to use CONS")
            self.scwlistdata=self.scw_data_nrt()
            return

        if self.scwversion=="any":
            if len(scw_cons)>0: # pick overlaps
                print("instructed to use ANY, and CONS is GOOD")
                self.scwlistdata=scw_cons
                return
            else:
                print("instructed to use ANY, and CONS is empty")
                self.scwlistdata=self.scw_data_nrt()
                return
        
        if self.scwversion=="latest":
            self.scwlistdata=self.scw_data_latest()



            

    def extract_from_index(self,scw_index,rep_base_prod=os.environ['INTEGRAL_DATA'],scwversion="001"):
        scx=SkyCoord(scw_index['RA_SCX'],scw_index['DEC_SCX'],unit="deg")

        target=SkyCoord(self.coordinates['RA'],self.coordinates['DEC'],unit="deg")

        m_avail=scw_index['SW_TYPE']=="POINTING"

        print(scw_index)

        m_coord=scx.separation(target).deg<self.coordinates['radius']

        print("for coordinates",target,"found",sum(m_coord & m_avail),"within",self.coordinates['radius'],"degrees")
        
        t1_ijd=float(converttime("UTC",self.timespan['T1'],"IJD"))
        t2_ijd=float(converttime("UTC",self.timespan['T2'],"IJD"))

        t1=scw_index['TSTART']
        t2=scw_index['TSTOP']
        m_time=(t1>t1_ijd) & (t2<t2_ijd)

        print("for time",t1_ijd,t2_ijd,"found",sum(m_coord & m_avail & m_time), "total range was",t1.min(),t2.max(),"(",(t2.max()-t1.min())/365,"years)")
            
        pre_selection=scw_index['SWID'][m_coord & m_avail & m_time]
        selection=[]

        if not self.randomize_pick:
            for scwid in pre_selection:
                evtsfn=rep_base_prod+"/scw/%s/%s.%s/isgri_events.fits.gz"%(scwid[:4],scwid,scwversion)
                print("searching for",evtsfn)
                if not os.path.exists(evtsfn):
                    print("skipping",scwid)
                    continue

                if not isfile(evtsfn) or not access(evtsfn, R_OK):
                    print("File {} doesn't exist or isn't readable".format(evtsfn))
                    continue

                selection.append(scwid)

                if self.max_pointings is not None and len(selection)>=self.max_pointings:
                    print("choosing only first",self.max_pointings)
                    break
        else:
            pick_size=min(self.max_pointings,len(pre_selection))
            print("choosing only random",pick_size)
            random.seed(0)
            pre_selection=sorted(random.sample(pre_selection,pick_size))
            
            for scwid in pre_selection:
                evtsfn=rep_base_prod+"/scw/%s/%s.%s/isgri_events.fits.gz"%(scwid[:4],scwid,scwversion)
                print("searching for",evtsfn)
                if not os.path.exists(evtsfn):
                    print("skipping",scwid)
                    continue

                if not isfile(evtsfn) or not access(evtsfn, R_OK):
                    print("File {} doesn't exist or isn't readable".format(evtsfn))
                    continue
                else:
                    print("File {} exist and readable".format(evtsfn))

                selection.append(scwid)


        print("selection:",selection)

        return [ddosa.ScWData(input_scwid=scwid+"."+scwversion) for scwid in selection]

class ReportScWList(ddosa.DataAnalysis):
    input_scwlist=TimeDirectionScWList

    cached=True

    def main(self):
        self.scwidlist=[scw.input_scwid.str() for scw in self.input_scwlist.scwlistdata]
