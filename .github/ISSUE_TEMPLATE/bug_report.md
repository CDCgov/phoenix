---
name: Bug report
about: Create a report to help us improve
title: "[BUG] - tile of issue"
labels: bug
assignees: ''

---

<!--- TO DO:
1.  Fill out the information below 
2. Delete the template and information not used 
3. Change the issue name  --->

**Describe the bug**
A clear and concise description of what feature is not working.

**Impact**
Please describe the impact this bug is causing.

**To Reproduce**
Steps to reproduce the behavior:
- What environment you are running the pipeline in (i.e. HPC, local, cloud etc).
- What version of the pipeline you are using
- The command that was run to cause the error
- If you used a custom config file please provide it
- The text of the error itself - screenshots are great for this!  
- Include the files that caused the error  
  - (if the file(s) contains something you don't want public open an issue with a brief description AND then email HAISeq@cdc.gov with the title QH ISSUE # so we can track the problem. DO NOT send PII!

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Logs**
If applicable, please attach logs to help describe your problem. These would be the `.command.err`, `.command.out`, `.command.sh` and/or the `.nextflow.log` files associated with the run/process that failed. The `.nextflow.log`  file is in directory you ran the pipeline from. The `.command.XXX` files are found in `<directory you ran pipeline from>/work/XX/xxxxxxxxxxxxxx` here the x and X are a random string of letters/numbers associated with the process. As the pipeline runs on CLI you will see the beginning of these strings to the left of the process that is running. 

**Additional context**
Add any other context about the problem here. If you have done any internet sleuthing already please link to any relevant posts on the topic.
