import urllib.request
import os
import time
from multiprocessing import Pool
import numpy as np
import argparse
import sys


def download_worker(url, timeout, header):
    try:
        header["user-agent"] = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.88 Safari/537.36"
        header["accept"] = "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9"
        this_request = urllib.request.Request(url, headers=header, method="GET")
        with urllib.request.urlopen(this_request) as f:
            content = f.read()
            print("finish")
            print(content)
            return content
    except Exception as e:
        print("Error(download_worker : {}): {}".format(url, e))
        try:
            error_code = e.args[0].split(" ", 1)[0]
            if error_code == "450":
                return  # no file
            elif error_code == "530":
                time.sleep(60)
                return download_worker(url, timeout, header)
            else:
                print(e.args)
                time.sleep(5)
                return download_worker(url, timeout, header)
        except:
            time.sleep(5)
            return download_worker(url, timeout, header)


def download(out_dir, info, worker, timeout, cache):
    filepath = "{}/{}".format(out_dir, info["filename"])
    task_number = np.ceil(info["size"] / cache).astype(np.int)
    pool = Pool(processes=worker)
    result_list = []
    if not os.path.exists(filepath):
        with open(filepath, "wb") as f:
            for ii in range(task_number):
                length = int(np.minimum((ii + 1) * cache, info["size"]) - ii * cache)
                contents = bytes("." * length, encoding="utf-8")
                f.write(contents)
                temp_size = os.path.getsize(filepath)
                done = int(50 * temp_size / info["size"])
                sys.stdout.write("\r[%s%s] %d%%" % ("█" * done, " " * (50 - done), 100 * temp_size / info["size"]))
                sys.stdout.flush()
        print()
    temp_size = os.path.getsize(filepath)
    if temp_size != info["size"]:
        with open(filepath, "wb") as f:
            write_size = (info["size"] - temp_size)
            for ii in range(np.ceil(write_size / cache).astype(np.int)):
                length = int(np.minimum((ii + 1) * cache, write_size) - ii * cache)
                contents = bytes("." * length, encoding="utf-8")
                f.write(contents)
                temp_size = os.path.getsize(filepath)
                done = int(50 * temp_size / info["size"])
                sys.stdout.write("\r[%s%s] %d%%" % ("█" * done, " " * (50 - done), 100 * temp_size / info["size"]))
                sys.stdout.flush()
        print()
    for ii in range(task_number):
        header = {"Range": "bytes={}-{}".format(ii * cache, np.minimum((ii + 1) * cache, info["size"]))}
        result_list.append(pool.apply_async(download_worker, args=(info["url"], timeout, header)))
    pool.close()
    result_index = list(range(len(result_list)))
    #  Same as .join()
    while len(list(filter(lambda x: not x.ready(), result_list))):
        ready_index = []
        for this_index, jj in enumerate(result_list):
            if jj.ready():
                ready_index.append(this_index)
        for this_index in ready_index[::-1]:
            ready_task = result_list.pop(this_index)
            ready_index = result_index.pop(this_index)
            with open(filepath, "ab") as f:
                f.seek(ready_index * cache)
                f.write(ready_task.get())
        temp_size = task_number - len(result_index)
        done = int(50 * temp_size / task_number)
        sys.stdout.write("\r[%s%s] %d%%" % ("█" * done, " " * (50 - done), 100 * temp_size / task_number))
        sys.stdout.flush()
    for ready_task, this_rindex in zip(result_list, result_index):
        with open(filepath, "ab") as f:
            f.seek(this_rindex * cache)
            f.write(ready_task.get())
    print("FINISH")


def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument("--url", required=True, type=str, help="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109762&format=file,https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124557&format=file")
    #parser.add_argument("--url", required=False, type=str, default="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109762&format=file,https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124557&format=file")
    parser.add_argument("--url", required=False, type=str, default="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3536499&format=file&file=GSM3536499%5FH1%5FK4me1%5FRep2%2Ebed%2Egz")
    #parser.add_argument("--out-dir", required=True, type=str, help="output folder")
    parser.add_argument("--out-dir", required=False, type=str, default="E:/monkey_single_cell/ChIP")
    parser.add_argument("--worker", required=False, type=int, default=5, help="")
    parser.add_argument("--timeout", required=False, type=int, default=30, help="")
    #parser.add_argument("--cache", required=False, type=int, default=512000, help="")
    parser.add_argument("--cache", required=False, type=int, default=512, help="")
    FLAGS = vars(parser.parse_args())
    url_list = FLAGS["url"].split(",")
    print("URL: {}".format(url_list))
    print("Output Dir: {}".format(FLAGS["out_dir"]))
    os.makedirs(FLAGS["out_dir"], exist_ok=True)
    download_list = []
    for this_url in url_list:
        this_request = urllib.request.Request(this_url)
        with urllib.request.urlopen(this_request) as f:
            content_info = {}
            for ii in [x.strip() for x in f.getheader("Content-Disposition").split(";")]:
                ii_fragments = ii.split("=")
                if len(ii_fragments) == 2:
                    content_info[ii_fragments[0]] = ii_fragments[1]
        assert "filename" in content_info.keys()
        assert "size" in content_info.keys()
        download_list.append({"filename": content_info["filename"].strip("\""),
                              "size": int(content_info["size"]),
                              "url": this_url})
    for ii in download_list:
        download(FLAGS["out_dir"], ii, FLAGS["worker"], FLAGS["timeout"], FLAGS["cache"])


if __name__ == "__main__":
    main()






