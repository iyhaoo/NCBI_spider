import urllib.request
import os
import time
import re
from ftplib import FTP
import json
import pandas as pd
import gzip
from multiprocessing import Pool
import threading

def grab_requset_key(request_data, timeoutLen=30, retry=0):
    try:
        request_type = re.search('[a-zA-Z]+?(?=[0-9])', request_data).group(0).lower()
        with urllib.request.urlopen("https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=" + request_data, timeout=timeoutLen) as urlOpen:
            response = urlOpen.read()
            search_result = re.search('(?<=RESULT={key:")[\w]+?(?=", mode:"run"})', str(response, encoding="utf-8"))
        if search_result:
            return request_type, search_result.group(0)
        else:
            if request_type == "gse":
                return request_type, "NA"
            else:
                print("Cannot find request key of {}, retry ...".format(request_data))
                return grab_requset_key(request_data, timeoutLen, retry - 1)
    except Exception as e:
        print("Error(grab_requset_key): {}".format(e))
        return grab_requset_key(request_data)

def request_summary(request_key, request_type, timeoutLen=30):
    try:
        if request_type == "srx":
            data = "stats=true&start=0&rows=1&fl=Assay_Type_s,AvgSpotLen_l,BioProject_s,BioSample_s,Center_Name_s,Consent_s,DATASTORE_filetype_ss,DATASTORE_provider_ss,Experiment_s,InsertSize_l,Instrument_s,LibraryLayout_s,LibrarySelection_s,LibrarySource_s,LoadDate_s,MBases_l,MBytes_l,Organism_s,Platform_s,ReleaseDate_s,Run_s,SRA_Sample_s,SRA_Study_s,Sample_Name_s,age_s,genotype_variation_s,source_name_s,strain_background_s,tissue_s&q=recordset:" + request_key
        else:
            data = "stats=true&start=0&rows=1&fl=Assay_Type_s,AvgSpotLen_l,BioProject_s,BioSampleModel_s,Center_Name_s,Consent_s,DATASTORE_provider_ss,InsertSize_l,LibraryLayout_s,LibrarySource_s,Organism_s,Platform_s,ReleaseDate_s,SRA_Study_s,age_s,biomaterial_provider_s,isolate_s,sex_s&q=recordset:" + request_key + "&stats.field=MBases_l&stats.field=MBytes_l"
        url = "https://trace.ncbi.nlm.nih.gov/Traces/study/proxy/run_selector.cgi?wt=json&indent=true&omitHeader=true&"
        req = urllib.request.Request(url=url, data=data.encode('utf-8'), method='POST')
        with urllib.request.urlopen(req, timeout=timeoutLen) as urlOpen:
            response = urlOpen.read()
            return_dict = json.loads(str(response, encoding="utf-8"))
            numFound = return_dict["response"]["numFound"]
            return_df = pd.DataFrame.from_dict(return_dict["response"]["docs"])
        return numFound, return_df
    except Exception as e:
        print("Error(request_summary): {}".format(e))
        return request_summary(request_key, request_type, timeoutLen)

def request_list_worker(req, timeoutLen):
    try:
        with urllib.request.urlopen(req, timeout=timeoutLen) as urlOpen:
            response = urlOpen.read()
            return_dict = json.loads(str(response, encoding="utf-8"))
        return pd.DataFrame.from_dict(return_dict["response"]["docs"])
    except Exception as e:
        print("Error(request_list_worker): {}".format(e))
        return request_list_worker(req, timeoutLen)

def request_list(request_key, numFound, timeoutLen=30):
    try:
        url = "https://trace.ncbi.nlm.nih.gov/Traces/study/proxy/run_selector.cgi?wt=json&indent=true&omitHeader=true&"
        nrows = 50
        return_df = pd.DataFrame()
        for ii in range(int(1 + numFound // nrows)):
            data = 'start={}&rows={}&fl=Run_s, BioSample_s, Sample_Name_s, DATASTORE_filetype_ss, AssemblyName_s, Experiment_s, Instrument_s, LibrarySelection_s, LoadDate_s, MBases_l, MBytes_l, SRA_Sample_s, cell_line_s, tissue_s&q=recordset:"{}"&sort=Run_s desc'.format(ii * nrows, nrows, request_key)
            req = urllib.request.Request(url=url, data=data.encode('utf-8'), method='POST')
            return_df = request_list_worker(req, timeoutLen) if return_df.empty else return_df.append(request_list_worker(req, timeoutLen))
            print("Extracted {}/{}".format(min(ii * nrows + nrows, numFound), numFound))
        return return_df
    except Exception as e:
        print("Error(request_list): {}".format(e))
        return request_list(request_key, numFound, timeoutLen)

def srx_request_worker(srx):
    request_type, request_key = grab_requset_key(srx)
    _, summary_df = request_summary(request_key, request_type)
    return summary_df


def soft_gz_request(soft_url, timeoutLen):
    try:
        with urllib.request.urlopen(soft_url, timeout=timeoutLen) as f:
            soft_gz = f.read()
        return str(gzip.decompress(soft_gz), encoding="utf-8")
    except:
        time.sleep(5)
        return soft_gz_request(soft_url, timeoutLen)

def get_srr_info_from_geo(dataset, maxThreadNum, timeoutLen):
    print("get_srr_info_from_geo")
    soft_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/{}nnn/{}/soft/{}_family.soft.gz".format(dataset[0:len(dataset) - 3], dataset, dataset)
    soft_gz = soft_gz_request(soft_url, timeoutLen)
    srx_list = [x.rsplit("=", 1)[1] for x in re.findall('(?<=!Sample_relation = SRA: ).+?SRX[0-9]+?(?=\n)', soft_gz)]
    p = Pool(maxThreadNum)
    result_list = []
    print("request srx info")
    for srx in srx_list:
        result = p.apply_async(srx_request_worker, kwds={"srx": srx})
        result_list.append(result)
    p.close()
    p.join()
    srx_all_df = pd.DataFrame()
    numFound = len(result_list)
    for ii, result in enumerate(result_list):
        srx_all_df = result.get() if srx_all_df.empty else srx_all_df.append(result.get())
        extracted_num = ii + 1
        if extracted_num % 50 == 0 or extracted_num == len(result_list):
            print("Extracted {}/{}".format(extracted_num, numFound))
    return srx_all_df

def just_make_summary(dataset, maxThreadNum=10, timeoutLen=30):
    request_type, request_key = grab_requset_key(dataset, timeoutLen)
    if request_key == "NA":
        srx_all_df = get_srr_info_from_geo(dataset, maxThreadNum, timeoutLen)
        out_summary = srx_all_df.to_csv(sep="\t", index=False)
    else:
        numFound, summary_df = request_summary(request_key, request_type, timeoutLen)
        list_df = request_list(request_key, numFound, timeoutLen)
        out_summary = "{}\n{}".format(summary_df.to_csv(sep="\t", index=False), list_df.to_csv(sep="\t", index=False))
    with open("{}/{}_summary.tsv".format(output_dir, dataset), "w") as f:
        f.write(out_summary)
    print("Extracted {} summary".format(dataset))

def checkpoint_make_or_load_srrUrl(dataset, tmp_dir, maxThreadNum=10, timeoutLen=30):
    srrUrl_file = "{}/{}_srrs.tsv".format(tmp_dir, dataset)
    srrUrl_list = []
    if os.path.exists(srrUrl_file):
        f = open(srrUrl_file)
        first_line = f.readline().strip()
        if first_line == "EOF":
            for line in f:
                srrUrl_list.append(line.strip())
            print("Loading infomation about {} ...".format(dataset))
            f.close()
            return srrUrl_list
        else:
            f.close()
    print("Searching infomation about {} ...".format(dataset))
    request_type, request_key = grab_requset_key(dataset, timeoutLen)
    if request_key == "NA":
        srx_all_df = get_srr_info_from_geo(dataset, maxThreadNum, timeoutLen)
        runs_list = srx_all_df["Run_s"].tolist()
        out_summary = srx_all_df.to_csv(sep="\t", index=False)
    else:
        numFound, summary_df = request_summary(request_key, request_type, timeoutLen)
        list_df = request_list(request_key, numFound, timeoutLen)
        runs_list = list_df["Run_s"].tolist()
        out_summary = "{}\n{}".format(summary_df.to_csv(sep="\t", index=False), list_df.to_csv(sep="\t", index=False))
    with open("{}/{}/{}_summary.tsv".format(output_dir, dataset, dataset), "w") as f:
        f.write(out_summary)
    with open(srrUrl_file, "w") as f:
        f.write("12345")  ###5 characters
        for SRRNO in runs_list:
            srrUrl = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{}/{}/{}.sra\n".format(
                SRRNO[:6], SRRNO, SRRNO)
            f.write(srrUrl)
            srrUrl_list.append(srrUrl.strip())
        f.seek(0)
        f.write("EOF\n")  ###5 characters
    return srrUrl_list

def main_threading(srr_name, run_info):
    try:
        if run_info["dataLen"]:
            print("Continue download {}: {}/{}".format(srr_name, run_info["finLen"], run_info["dataLen"]))
        else:
            print("Begin download {} ...".format(srr_name))
        #print(run_info)
        srrUrl_split = run_info["srrUrl"].split("/", 3)
        ftpServer = srrUrl_split[2]
        remotePath = srrUrl_split[3].strip()
        timeoutLen = run_info["timeoutLen"]
        ftp = FTP(host=ftpServer, timeout=timeoutLen)
        ftp.login()
        run_info["cache"] = ""
        if not run_info["dataLen"]:
            def ensure_data_length(self):
                regex = re.compile('[a-zA-Z]+?[ ]+?[0-9]+?(?=[ ]+?[a-zA-Z]+?)')
                run_info["dataLen"] = int(str(regex.findall(str(self))[0]).split(" ")[-1])
            ftp.retrlines('LIST ' + remotePath, callback=ensure_data_length)
        def ftp_callback(self):
            if len(run_info["cache"]) == 0:
                run_info["cache"] = self
            else:
                run_info["cache"] += self
            if len(run_info["cache"]) >= run_info["cacheLen"] or run_info["finLen"] + len(run_info["cache"]) >= run_info["dataLen"]:
                file_path = "{}/{}".format(run_info["save_dir"], srr_name)
                f = open(file_path, "wb") if run_info["finLen"] == 0 else open(file_path, "ab")
                f.write(run_info["cache"])
                f.close()
                run_info["finLen"] += len(run_info["cache"])
                use_time_once = max(time.time() - run_info["last_time"], 1e-8)
                print("wrinting {}\t{}/{} ({:.2%})\tSpeed: {}k/s"
                      .format(file_path, run_info["finLen"], run_info["dataLen"],
                              run_info["finLen"] / run_info["dataLen"],
                              round(len(run_info["cache"]) / (1000 * use_time_once), 2)))
                run_info["cache"] = ""
                if os.path.getsize(file_path) != run_info["finLen"]:
                    os.remove(file_path)
                    run_info["finLen"] = 0
                    raise NameError(file_path, "Error: Inconsistent file and log")

                with open("{}/{}.log".format(run_info["log_dir"], srr_name.split(".", 1)[0]), "w") as log_f:
                    log_f.write("{}\t{}\t{}\n".format(srr_name, run_info["dataLen"], run_info["finLen"]))
                run_info["last_time"] = time.time()
        ftp.retrbinary(cmd='RETR ' + remotePath, callback=ftp_callback, rest=run_info["finLen"])
        ftp.quit()
        print("Finish: {} (Use {} seconds)".format(srr_name, round(time.time() - run_info["first_time"], 2)))
    except Exception as e:
        print("Error(main_threading of {}): {}".format(srr_name, e))
        try:
            error_code = e.args[0].split(" ", 1)[0]
            if error_code == "450":
                return  # no file
            elif error_code == "530":
                time.sleep(60)
                return main_threading(srr_name, run_info)
            else:
                print(e.args)
                time.sleep(5)
                return main_threading(srr_name, run_info)
        except:
            time.sleep(5)
            return main_threading(srr_name, run_info)

class ftpFileDownloadThread(threading.Thread):
    def __init__(self, run_item):
        threading.Thread.__init__(self)
        self.run_item = run_item
    def run(self):
        srr_name, run_info = self.run_item
        run_info["last_time"] = time.time()
        run_info["first_time"] = run_info["last_time"]
        main_threading(srr_name, run_info)

def ftpFileDownload(srrUrl_dict_items, output_dir, maxThreadNum, timeoutLen=30, cacheLen=512000):
    log_dir = "{}/tmp/downloadLog".format(output_dir)
    dataset, srrUrl_list = srrUrl_dict_items
    run_info_dict = {}
    for srrUrl in srrUrl_list:
        srr_name = srrUrl.rsplit("/", 1)[1]
        run_info_dict[srr_name] = {"dataLen": None,
                                   "finLen": 0,
                                   "last_time": time.time(),
                                   "srrUrl": srrUrl,
                                   "save_dir": "{}/{}".format(output_dir, dataset),
                                   "log_dir": log_dir,
                                   "cacheLen": cacheLen,
                                   "timeoutLen": timeoutLen}
    if os.path.isdir(log_dir):
        for log_file in os.listdir(log_dir):
            with open("{}/{}".format(log_dir, log_file)) as f:
                downloadState = f.read().replace("\n", "").split("\t")
                if downloadState[0] in run_info_dict:
                    if int(downloadState[1]) == int(downloadState[2]):
                        del run_info_dict[downloadState[0]]
                    else:
                        run_info_dict[downloadState[0]]["dataLen"] = int(downloadState[1])
                        run_info_dict[downloadState[0]]["finLen"] = int(downloadState[2])

    else:
        os.mkdir(log_dir)
    for run_item in run_info_dict.items():
        thisThread = ftpFileDownloadThread(run_item)
        while threading.activeCount() > maxThreadNum:
            time.sleep(1)
        thisThread.start()









#
#
#
#
#
# download settings#######

output_dir = "E:/single_cell_data/ncbi_spider_20181022"
allDatasets = ["GSE96772", "GSE110513", "GSE75478"]
#just_summary = True
just_summary = False
#
#
#
#
#
#  settings#######
maxThreadNum = 10
timeoutLen = 30
cacheLen = 512000
#  settings_end###
#
#
#
#
#
if __name__ == '__main__':
    if just_summary:
        os.makedirs(output_dir, exist_ok=True)
        for dataset in allDatasets:
            just_make_summary(dataset, maxThreadNum, timeoutLen)
    else:
        tmp_dir = output_dir + "/tmp"
        os.makedirs(tmp_dir, exist_ok=True)
        print("Dest_dir: {}".format(output_dir))
        srrUrl_dict = {}
        for index, dataset in enumerate(allDatasets):
            output_dataset_dir = "{}/{}".format(output_dir, dataset)
            os.makedirs(output_dataset_dir, exist_ok=True)
            srrUrl_dict[dataset] = checkpoint_make_or_load_srrUrl(dataset, tmp_dir, maxThreadNum, timeoutLen)
        file_thread_dict = {}
        for srrUrl_dict_items in srrUrl_dict.items():
            ftpFileDownload(srrUrl_dict_items, output_dir, maxThreadNum, timeoutLen, cacheLen)

