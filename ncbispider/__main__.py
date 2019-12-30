import urllib.request
import os
import time
import re
from ftplib import FTP
import pandas as pd
import gzip
from multiprocessing import Pool
import threading
import argparse


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
    soft_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/{}nnn/{}/soft/{}_family.soft.gz".format(dataset[0:len(dataset) - 3], dataset, dataset)
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

def just_make_summary(dataset, output_dir, maxThreadNum, timeoutLen):
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

def checkpoint_make_or_load_srrUrl(dataset, tmp_dir, output_dir, maxThreadNum, timeoutLen):
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

def main_threading(file_name, run_info):
    try:
        if run_info["dataLen"]:
            print("Continue download {}: {}/{}".format(file_name, run_info["finLen"], run_info["dataLen"]))
        else:
            print("Begin download {} ...".format(file_name))
        url_fragments = run_info["srrUrl"].split("/", 3)
        ftpServer = url_fragments[2]
        remotePath = url_fragments[3].strip()
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
                file_path = "{}/{}".format(run_info["save_dir"], file_name)
                f = open(file_path, "wb") if run_info["finLen"] == 0 else open(file_path, "ab")
                f.write(run_info["cache"])
                f.close()
                run_info["finLen"] += len(run_info["cache"])
                use_time_once = max(time.time() - run_info["last_time"], 1e-8)
                print("wrinting {}\t{}/{} ({:.2%})\tSpeed: {}k/s".format(file_path,
                                                                         run_info["finLen"],
                                                                         run_info["dataLen"],
                                                                         run_info["finLen"] / run_info["dataLen"],
                                                                         round(len(run_info["cache"]) / (1000 * use_time_once), 2)))
                run_info["cache"] = ""
                if os.path.getsize(file_path) != run_info["finLen"]:
                    os.remove(file_path)
                    run_info["finLen"] = 0
                    raise NameError(file_path, "Error: Inconsistent file and log")
                with open("{}/{}.log".format(run_info["log_dir"], file_name.split(".", 1)[0]), "w") as log_f:
                    log_f.write("{}\t{}\t{}\n".format(file_name, run_info["dataLen"], run_info["finLen"]))
                run_info["last_time"] = time.time()
        ftp.retrbinary(cmd='RETR ' + remotePath, callback=ftp_callback, rest=run_info["finLen"])
        ftp.quit()
        print("Finish: {} (Use {} seconds)".format(file_name, round(time.time() - run_info["first_time"], 2)))
    except Exception as e:
        print("Error(main_threading of {}): {}".format(file_name, e))
        try:
            error_code = e.args[0].split(" ", 1)[0]
            if error_code == "450":
                return  # no file
            elif error_code == "530":
                time.sleep(60)
                return main_threading(file_name, run_info)
            else:
                print(e.args)
                time.sleep(5)
                return main_threading(file_name, run_info)
        except:
            time.sleep(5)
            return main_threading(file_name, run_info)

class ftpFileDownloadThread(threading.Thread):
    def __init__(self, run_item):
        threading.Thread.__init__(self)
        self.run_item = run_item
    def run(self):
        file_name, run_info = self.run_item
        run_info["last_time"] = time.time()
        run_info["first_time"] = run_info["last_time"]
        main_threading(file_name, run_info)

def ftpFileDownload(srrUrl_dict_items, output_dir, maxThreadNum, timeoutLen, cacheLen):
    log_dir = "{}/tmp/downloadLog".format(output_dir)
    dataset, srrUrl_list = srrUrl_dict_items
    run_info_dict = {}
    for srrUrl in srrUrl_list:
        file_name = srrUrl.rsplit("/", 1)[1]
        run_info_dict[file_name] = {"dataLen": None,
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=str, help="GSE96772,GSE110513,GSE75478")
    parser.add_argument("--out-dir", required=True, type=str, help="output folder")
    parser.add_argument("--maxThreadNum", required=False, type=int, default=10, help="")
    parser.add_argument("--timeoutLen", required=False, type=int, default=30, help="")
    parser.add_argument("--cacheLen", required=False, type=int, default=512000, help="")
    FLAGS = vars(parser.parse_args())
    dataset_list = FLAGS["dataset"].split(",")
    print("Dataset: {}".format(dataset_list))
    print("Output Dir: {}".format(FLAGS["out_dir"]))
    tmp_dir = FLAGS["out_dir"] + "/tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    srrUrl_dict = {}
    for index, dataset in enumerate(dataset_list):
        output_dataset_dir = "{}/{}".format(FLAGS["out_dir"], dataset)
        os.makedirs(output_dataset_dir, exist_ok=True)
        srrUrl_dict[dataset] = checkpoint_make_or_load_srrUrl(dataset, tmp_dir, FLAGS["out_dir"], FLAGS["maxThreadNum"], FLAGS["timeoutLen"])
    for srrUrl_dict_items in srrUrl_dict.items():
        ftpFileDownload(srrUrl_dict_items, FLAGS["out_dir"], FLAGS["maxThreadNum"], FLAGS["timeoutLen"], FLAGS["cacheLen"])


if __name__ == "__main__":
    main()






