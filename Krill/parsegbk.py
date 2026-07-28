#!/usr/bin/env python
# Copyright (C) 2015,2016 Mohammad Alanjary
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
# This file is part of ARTS
# ARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
# License: You should have received a copy of the GNU General Public License v3 with ARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import os
from Bio import SeqIO
import setlog
import unicodedata
print("\n \n ################# ARTS ################################## ARTS #################")
print("\n \n este arquivo esta sendo executado: parsegbk.py")
print("\n \n ################# ARTS ################################## ARTS #################")
global log
log = setlog.init(toconsole=True)

def getheader(seq_record,recnum,userecnum=False):
    """Get record header information"""
    if len(seq_record.id) and not userecnum:
        seqtitle = seq_record.id
    elif len(seq_record.name) and not userecnum:
        seqtitle = seq_record.name
    else:
        log.info("record has no ID using reccord number")
        seqtitle = "scaffold_"+str(recnum)

    if len(seq_record.description) > 0:
        seqdesc = seq_record.description
    elif "source" in seq_record.annotations.keys():
        seqdesc = seq_record.annotations["source"].replace(".","").replace(" ","_")

    else:
        seqdesc = ""
        log.info("Genbank does not have record description")
    return seqtitle, seqdesc


def appendheader(seqtitle, desc, lnum, qual, loc, recnum):
    """Add details to fasta header"""
    gi = "lcl|" + str(lnum)
    gdesc = "-"
    gdesc2 = "-"
    if "db_xref" in qual.keys():
        gi = qual["db_xref"][0].replace("|", "-").lower()
        if ":" in gi:
            gi = gi.replace(":", "|")
    if "product" in qual.keys():
        gdesc = qual["product"][0]
    if "gene" in qual.keys():
        gdesc2 = qual["gene"][0]
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("gi line 90: \n", gi)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("seqtitle line 90: \n", seqtitle)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("desc line 90: \n", desc)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("gdesc2 line 90: \n", gdesc2)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("gdesc line 90: \n", gdesc)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("loc[0] line 90: \n", loc[0])
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("loc[1] line 90: \n", loc[1])
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("loc[2] line 90: \n", loc[2])
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("recnum line 90: \n", recnum)
        print("\n ############## ARTS parsegbk.py ########### \n")
        print("lnum line 90: \n", lnum)
    return "%s|%s %s|%s|%s|loc|%s %s %s # ID=%s_%s"%\
           (gi,seqtitle,desc,gdesc2,gdesc,loc[0],loc[1],loc[2],recnum,lnum)


    #return gi + "|" + seqtitle + " " + desc + "|" + gdesc2 + "|" + gdesc + "|loc|" + str(loc[0]) + " " + str(loc[1] + 1) + " " + str(loc[2]) + " # ID=1_" + str(lnum) + ";"


def convertgenes(filename, outdir="./",rename=False,usetrans=False,plasmid=False,userecnum=False,clust=False,cutoff=10):
    """Parse all gbk records and output nuc and prot sequences for each CDS in multi-fasta format"""
    log.info("Starting %s..."%filename)
    fpath, fname = os.path.split(filename)
    fname, ext = os.path.splitext(fname)
    reclist = SeqIO.parse(filename, "genbank")
    genus = "unknown"
    orgname = "unknown"

    if not outdir.endswith("/"):
        outdir+="/"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    lnum = 1
    clusters = 0
    cdscount = 0
    recnum = 0
    features=["CDS","rRNA"]
    clusters=[]
    with open(outdir + fname + ".clust.tsv", "w") as clust_handle, open(outdir + fname + ".fna", "w") as nuc_handle:
        clust_handle.write("#cluster_number\tproduct\tscaffold/record_number\tstart\tend\n")
        for seq_record in reclist:
            recnum += 1
            seqtitle, seqdesc = getheader(seq_record,recnum,userecnum)
            print("\n ############## ARTS parsegbk.py ########### \n")
            print("seqtitle line 119: \n", seqtitle)
            print("\n ############## ARTS parsegbk.py ########### \n")
            print("seqdesc line 119: \n", seqdesc)
            bpcount = 0
            bpend = 1
            for seq_feature in seq_record.features:
                plasmidTitle=""
                if "source" in seq_feature.type.lower() and "plasmid" in seq_feature.qualifiers.keys():
                    plasmidTitle="_PLASMID_"+str(seq_feature.qualifiers["plasmid"])
                if seq_feature.type in features:
                    outseq = seq_feature.extract(seq_record.seq)
                    if len(outseq) > cutoff:
                        bpend = int(seq_feature.location.end)
                        bpstart = int(seq_feature.location.start)
                        #Fix for CDS joins that span multiple strands (defaults to 1)
                        if seq_feature.location.strand:
                            bpstrand = int(seq_feature.location.strand)
                        else:
                            bpstrand = 1
                        bpcount += bpend - bpstart
                        seqdetails = appendheader(seqtitle, seqdesc, lnum, seq_feature.qualifiers,
                                                  [bpstart, bpend, bpstrand], recnum)
                        seqdetails+=plasmidTitle
                        seqdetails+="|Type="+seq_feature.type
                        if "transl_table" in seq_feature.qualifiers:
                            seqdetails += "|transl_table="+seq_feature.qualifiers["transl_table"][0]

                        #Remove non-unicode chars
                        #seqdetails = "".join([x if ord(x) < 128 else '-' for x in str(seqdetails)])
                        seqdetails = str(seqdetails).decode('ascii','ignore')
                        nuc_handle.write(">%s\n%s\n" % (seqdetails, outseq))
                        lnum += 1
                    if seq_feature.type == "CDS":
                        cdscount += 1

                #get cluster if present
                if "cluster" in seq_feature.type.lower():
                    note = seq_feature.qualifiers.get("note",[])
                    product = seq_feature.qualifiers.get("product",[])
                    bpend = int(seq_feature.location.end)
                    bpstart = int(seq_feature.location.start)
                    for x in note:
                        if "Cluster number:" in x:
                            clustnum = x.split(":")[-1].strip()
                            clusters.append([clustnum,str(",".join(product)),seqtitle,bpstart,bpend])
                            clust_handle.write("%s\t%s\t%s\t%s\t%s\n"%tuple(clusters[-1]))

                    # if not clustnum:
                    #     clustnum = clusters

            log.info("Record #%s CDS bp coverage: %s%%"%(recnum,(bpcount * 100 / bpend)))

            if "organism" in seq_record.annotations:
                orgname=seq_record.annotations["organism"]
                genus = orgname.replace("_"," ").split()[0]
            elif "source" in seq_record.annotations:
                orgname=seq_record.annotations["source"]
                genus = orgname.replace("_"," ").split()[0]
            else:
                orgname=""
                genus=""
    if rename:
        os.rename(outdir + fname + ".fna", outdir + rename + ".fna")
        os.rename(outdir + fname + ".clust.tsv", outdir + rename + ".clust.tsv")
        outfile=rename
    else:
        outfile=fname
    #remove any bad chars
    orgname = ''.join(ch for ch in orgname.replace(".","").replace(",","").replace("(","_").replace(")","_").replace("-","_").replace(" ","_") if ch.isalnum() or "_")
    genus = ''.join(ch for ch in genus.replace(".","").replace(",","").replace("(","_").replace(")","_").replace("-","_").replace(" ","_") if ch.isalnum() or "_")
    if not clust:
        os.remove(outdir + outfile + ".clust.tsv")
    if len(orgname) < 3:
        orgname = "unknown"
        genus = "unknown"
        log.warning("Could not find organism name or genus... setting to 'unknown'. Add SOURCE / ORGANISM sections to genbank to resolve")

    log.info("# Finished %s # Records: %s; CDS features: %s; Clusters: %s"%(filename,recnum,cdscount,len(clusters)))
    return orgname,genus,clusters