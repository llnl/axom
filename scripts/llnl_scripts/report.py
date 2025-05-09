#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: report.py

 description: 
  Emails an HTML report for jobs and test results

"""

import datetime
import operator
import os
import re
import shutil
import sys

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from argparse import ArgumentParser
from smtplib import SMTP
from email.mime.text import MIMEText
from dateutil.parser import parse
from llnl_lc_build_tools import *


__cellBlackOnGreen = "    <td bgcolor=\"#00C000\"><font color=\"#000000\">{0}</font></td>\n"
__cellWhiteOnRed   = "    <td bgcolor=\"#E10000\"><font color=\"#FFFFFF\">{0}</font></td>\n"
__cellNormal       = "    <td>{0}</td>\n"

__contentHeader = "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></head><body>\n"
__contentFooter = "</body>\n</html>\n"
__detailLinkIndex = 0


class JobInfo(object):
    def __init__(self):
        self.name = ""
        self.sys_type = ""
        self.archive_dir = ""
        self.success = False
        self.datetime = ""
        self.specInfos = {}


class SpecInfo(object):
    def __init__(self):
        self.name = ""
        self.success = False
        self.passed = []
        self.failed = []

def parse_args():
    "Parses args from command line"
    parser = ArgumentParser()
    
    parser.add_argument("--clean",
                        dest="clean",
                        action="store_true",
                        default=False,
                        help="!!DESTRUCTIVE!! This option cleans old archived jobs (leaves 20 jobs)")

    parser.add_argument("-e", "--email",
                        type=str,
                        dest="email",
                        default="axom-dev@llnl.gov",
                        help="Email address to send report (Defaults to 'axom-dev@llnl.gov')")

    parser.add_argument("--html",
                        type=str,
                        dest="html",
                        default="status.html",
                        help="File name for saving generated html status file (Defaults to 'status.html')")

    ###############
    # parse args
    ###############
    opts = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using argparse
    opts = vars(opts)
    return opts


def main():
    opts = parse_args()

    # Email info
    sender = "white238@llnl.gov"
    receiver = opts["email"]
    emailServer = "nospam.llnl.gov"
    if on_rz():
        emailSubject = "[AXOM] RZ Status"
    else:
        emailSubject = "[AXOM] CZ Status"

    archive_dir = get_archive_base_dir()

    if opts["clean"]:
        cleanOldArchives(archive_dir)

    #Generate build email and send
    print("Reading archived job information...")
    basicJobInfos, srcJobInfos, tplJobInfos = generateJobInfos(archive_dir)

    print("Generating email content...")
    emailContent = generateEmailContent(basicJobInfos, srcJobInfos, tplJobInfos)
    
    print("Saving html file '{}'".format( opts["html"] ))
    with open(opts["html"], 'w') as f:
        f.write(emailContent)

    print("Sending email to {}...".format(opts["email"]))
    return sendEmail(emailContent, emailSubject, sender, receiver, emailServer)


def getAllDirectoryNames(path):
    return [n for n in os.listdir(path) if os.path.isdir(pjoin(path, n))]


def cleanOldArchives(archive_dir):
    print("Deleting old archive directories...")

    # Remove only the individual jobs not any of the directory structure
    # Directory structure = <archive base>/<sys_type>/<job name>/<datetime>
    sys_types = getAllDirectoryNames(archive_dir)
    for sys_type in sys_types:
        sys_type_dir = pjoin(archive_dir, sys_type)
        jobNames = getAllDirectoryNames(sys_type_dir)
        for jobName in jobNames:
            jobName_dir = pjoin(sys_type_dir, jobName)
            datetimes = getAllDirectoryNames(jobName_dir)

            datetimes.sort()
            datetimes.reverse()
            del datetimes[:20]
            for datetime in datetimes:
                datetime_dir = pjoin(jobName_dir, datetime)
                shutil.rmtree(datetime_dir)

    print("Done deleting.")


def determineSuccessState(path):
    success = False
    if os.path.exists(pjoin(path, "success.json")):
        success = True
    elif os.path.exists(pjoin(path, "failed.json")):
        success = False
    else:
        success = None
    return success


def splitJobName(name):
    # Job Names should follow the following naming scheme "MachineName_-_JobName"
    if not "_-_" in name:
        machine = ""
        jobName = name
    else:
        splitName = name.split('_-_', 1)
        machine = splitName[0]
        jobName = splitName[1]
    return machine, jobName


def generateJobInfos(archive_dir):
    basicJobInfos = []
    srcJobInfos = []
    tplJobInfos = []
    sys_types = getAllDirectoryNames(archive_dir)

    # Directory structure = <archive base>/<sys_type>/<job name>/<datetime>/<spec>/[Logs]
    for sys_type in sys_types:
        sys_type_dir = pjoin(archive_dir, sys_type)
        jobNames = getAllDirectoryNames(sys_type_dir)
        for jobName in jobNames:
            # We only want the newest
            jobName_dir = pjoin(sys_type_dir, jobName)
            datetimes = getAllDirectoryNames(jobName_dir)
            datetimes.sort()
            datetime = datetimes[-1]
            datetime_dir = pjoin(jobName_dir, datetime)

            # Populate job information
            currJobInfo = JobInfo()
            currJobInfo.name = jobName
            currJobInfo.sys_type = sys_type
            currJobInfo.datetime = datetime
            currJobInfo.archive_dir = datetime_dir
            currJobInfo.success = determineSuccessState(datetime_dir)

            # If any specs add their information
            specs = getAllDirectoryNames(datetime_dir)
            for spec in specs:
                currSpecInfo = SpecInfo()
                currSpecInfo.name = spec
                spec_dir = pjoin(datetime_dir, spec)
                currSpecInfo.success = determineSuccessState(spec_dir)

                populateTests(currSpecInfo, spec_dir)
                if currJobInfo.specInfos.has_key(currSpecInfo.name):
                    print("Warning: duplicate spec ({0}) found in job: {1}".format(spec, currJobInfo.name))
                currJobInfo.specInfos[spec] = currSpecInfo

                currJobInfo.success = True
                for specName in currJobInfo.specInfos.keys():
                    # if any tests have failed then the whole job fails
                    currSpecInfo = currJobInfo.specInfos[specName]
                    if len(currSpecInfo.failed) > 0 or currSpecInfo.success is None:
                        currJobInfo.success = False
                        break
            
            # Add job to corresponding list
            if len(specs) == 0:
                # Basic jobs should contain no folders
                # Used for general purpose jobs
                basicJobInfos.append(currJobInfo)
            else:
                # TPL jobs will have spec dirs and an info.json at the root
                if os.path.exists(pjoin(datetime_dir, "info.json")):
                    tplJobInfos.append(currJobInfo)
                else:
                    srcJobInfos.append(currJobInfo)

    return basicJobInfos, srcJobInfos, tplJobInfos


def populateTests(specInfo, path):
    test_xml_path = pjoin(path, 'Test.xml')
    if not os.path.exists(test_xml_path):
        return

    #xml_root = ET.parse(test_xml_path).getroot()
    tree = ET.parse(test_xml_path)

    for test_elem in tree.iter(tag="Test"):
        # Skip the list of test names @ Site -> Testing -> TestList -> Test
        if not test_elem.attrib.has_key("Status"):
            continue

        name = ""
        for child_elem in test_elem:
            if child_elem.tag == "Name":
                name = child_elem.text
                break
        if name == "":
            print("Error: {0}: Unknown to find test name".format(test_xml_path))

        if test_elem.attrib["Status"] == "passed":
            specInfo.passed.append(name)
        elif test_elem.attrib["Status"] == "failed":
            specInfo.failed.append(name)
        else:
            print("Error: {0}: Unknown test status ({1})".format(test_xml_path, test_elem.attrib["Status"]))


def generateEmailContent(basicJobInfos, srcJobInfos, tplJobInfos):
    html = __contentHeader

    # Add Basic jobs to the email
    if len(basicJobInfos) > 0:
        # Header for all basic jobs
        html += "<center><font size=\"5\"> </font></center><br>\n"
        html += "<center><b><font size=\"5\">Basic Jobs</font></b></center>\n"
        html += "<center><table border=\"1\">\n"

        # Start table for basic jobs
        html += "<tr>\n"
        html += "    <th bgcolor=\"#C0C0C0\">Machine</th>\n"
        html += "    <th bgcolor=\"#C0C0C0\">Job</th>\n"
        html += "    <th bgcolor=\"#C0C0C0\">Date</th>\n"
        html += "    <th bgcolor=\"#C0C0C0\">Success</th>\n"
        html += "</tr>\n\n"

        basicJobInfos.sort(key=operator.attrgetter("name"))
        for currJobInfo in basicJobInfos:
            machine, jobName = splitJobName(currJobInfo.name)
            # Row for each job
            html += "<tr>\n"
            html += "    <td>" + machine + "</td>\n"
            html += "    <td>" + jobName + "</td>\n"
            html += getDateCellHtml(currJobInfo.datetime)
            html += getSuccessFormat(currJobInfo.success,
                                     __cellBlackOnGreen, __cellWhiteOnRed).format(str(currJobInfo.success))
            html += "</tr>\n\n"

        # Close table for basic jobs
        html += "</table></center>\n\n"

    html += getHTMLforJobInfos(srcJobInfos, False)
    html += getHTMLforJobInfos(tplJobInfos, True)

    # Add date time to the end of the email to make listserv not reject duplicate emails
    html += "<br/>Generated at " + get_timestamp()

    html += __contentFooter
    return html


def getHTMLforJobInfos(jobInfos, isTPLJob):
    if len(jobInfos) == 0:
        return ""

    html = ""

    # Header for all src jobs
    html += "<center><font size=\"5\"> </font></center><br>\n"
    if (isTPLJob):
        html += "<center><b><font size=\"5\">TPL Jobs</font></b></center>\n"
    else:    
        html += "<center><b><font size=\"5\">Source Jobs</font></b></center>\n"

    # Get and sort all sys types
    sys_types = []
    for currJobInfo in jobInfos:
        sys_types.append(currJobInfo.sys_type)
    sys_types = list(set(sys_types))
    sys_types.sort()

    # 2 = number of table cells that are not specs (Name, Date)
    baseCellWidth = 2

    # Pre calculate the table width
    tableWidth = 0
    for currJobInfo in jobInfos:
        specNames = currJobInfo.specInfos.keys()
        tableWidth = max(tableWidth, baseCellWidth + len(specNames))

    # Write table headers
    html += "<table border=\"1\">\n"
    html += "<tr>\n"
    html += "    <th bgcolor=\"#C0C0C0\">Machine</th>\n"
    html += "    <th bgcolor=\"#C0C0C0\">Date</th>\n"
    for count in range(0,tableWidth-baseCellWidth):
        html += "    <th bgcolor=\"#C0C0C0\"></th>\n"
    html += "</tr>\n\n"

    jobInfos.sort(key=operator.attrgetter("name"))
    for currJobInfo in jobInfos:
        # Row for each job
        html += "<tr>\n"
        html += "    <td>" + splitJobName(currJobInfo.name)[0] + "</td>\n"
        html += getDateCellHtml(currJobInfo.datetime)

        failedTestsString = ""
        # Write cell for each spec
        for specName in currJobInfo.specInfos.keys():
            currSpecInfo = currJobInfo.specInfos[specName]
            totalNumberOfTests = len(currSpecInfo.passed) + len(currSpecInfo.failed)
            cellString = "<pre>{0}\n{1}/{2}</pre>".format(specName.split('~')[0],
                                              str(len(currSpecInfo.passed)),
                                              str(totalNumberOfTests))
            html += getSuccessFormat(currSpecInfo.success,
                                     __cellBlackOnGreen, __cellWhiteOnRed).format(cellString)

            # Build up a list of failed tests to be added to the html table at the end
            if len(currSpecInfo.failed) > 0:
                if failedTestsString != "":
                    failedTestsString += "<br />"
                failedTestsString += "[{0}]".format(currSpecInfo.name)
                for name in currSpecInfo.failed:
                    if failedTestsString[-1] != "]":
                        failedTestsString += ",\n"
                    failedTestsString += " " + name

        # Pad table size if necessary
        widthDifference = tableWidth - baseCellWidth - len(currJobInfo.specInfos.keys())
        if widthDifference != 0:
            html += "    <td colspan=\"{0}\"></td>\n".format(str(widthDifference))
        html += "</tr>\n"

        if failedTestsString != "":
            html += "<tr><td colspan=\"{0}\">{1}</td></tr>\n\n".format(str(tableWidth), failedTestsString)
        html += "<tr><td colspan=\"{0}\">{1}</td></tr>\n\n".format(str(tableWidth),currJobInfo.archive_dir)

    # Close table
    html += "</table>\n"

    return html


def getDateCellHtml(datetimeString):
    readableDatetime = convertDatetimeToReadable(datetimeString)
    # Color cell if datetime is older than a week
    isOld = isDatetimeOld(datetimeString)
    html = getSuccessFormat(isOld, __cellNormal, __cellWhiteOnRed).format(readableDatetime)
    return html


def convertDatetimeToReadable(datetimestring):
    dt = datetime.datetime.strptime(datetimestring, "%Y_%m_%d_%H_%M_%S")
    dateString = dt.strftime("%a %m/%d/%y %H:%M")
    return dateString


def getSuccessFormat(success, successFormat, failureFormat):
    if success:
        return successFormat
    return failureFormat


def sendEmail(content, subject, sender, receiver, emailServer):
    try:
        msg = MIMEText(content, "html")
        msg["Subject"] = subject
        msg["From"] = sender
        msg["To"] = receiver
        msg["Content-Type"] = "text/html"

        conn = SMTP(emailServer)
        conn.set_debuglevel(False)
        try:
            conn.sendmail(sender, receiver, msg.as_string())
        finally:
            conn.close()
    except Exception as e:
        print("Failed to send email:\n {0}".format(str(e)))
        return False
    print("Sent.")

    return True


def isDatetimeOld(datetimestring):
    dt = datetime.datetime.strptime(datetimestring, "%Y_%m_%d_%H_%M_%S")
    isOld = False
    if datetime.timedelta(days=7) >= datetime.datetime.now()-dt:
        isOld = True
    return isOld


if __name__ == "__main__":
    if main():
        sys.exit(0)
    sys.exit(1)

