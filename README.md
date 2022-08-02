# :fire::bird::fire: PHoeNIx

<!-- [![GitHub Downloads](https://img.shields.io/github/downloads/CDCgov/phoenix/total.svg?style=social&logo=github&label=Download)](https://github.com/CDCgov/phoenix/releases) -->
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

For full documentation on the pipeline see the [Wiki](https://github.com/cdcent/phoenix/wiki), but quick start instructions are provided below if you are feeling brave. 

# Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`). 

   There are several options for install if you do not already have it on your system:

   * Install into conda environment, which will require a version of Anaconda to be installed on your system.

       ```console
       mamba create -n nextflow -c bioconda nextflow=21.10.6  
       ```

      <!---```console
       mamba create -n nextflow -c bioconda -c conda-forge nf-core=2.2 nextflow=21.10.6 git=2.35.0 openjdk=8.0.312 graphviz
       ```--->

   * If you prefer a to use `curl` or `wget` for install see the [Nextflow Documentaiton](https://www.nextflow.io/docs/latest/getstarted.html) 

2. Install [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity >=3.8.0`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility. 

3. Email HAISeq@cdc.gov, with the subject line "krakenDB invite request" to request access to the sharefile link and provide the email address to send invite to.

4. Download the `hash.k2d`, `opts.k2d`, and `taxo.k2d` files needed for the kraken2 subworkflow of PHoeNIx from the CDC sharefile link. You **CANNOT** use a different krakenDB for this as it needs to match the `ktax_map.k2` file that is included in the pipeline. At this time this is not downloadable via command line. Once downloaded the folder containing these files is passed to PHoeNIx via the `--kraken2db`.
    
     > If you ran the v1.0.0-dev version of the pipeline and already downloaded the `hash.k2d` file there is no need to redownload it. The `opts.k2d`, and `taxo.k2d` are found in the `phoenix/assets/databases` folder. You can just copy these over into a new folder to pass to the PHoeNIx.  

5. (optional) If you installed nextflow via a conda environment activate the nextflow environment with:  

   ```console
   conda activate nextflow
   ```

6. Run PHoeNIx on a test sample loaded with the package with a single command:

    ```console
    nextflow run cdcgov/phoenix -r v1.0.0 -profile <singularity/docker/custom>,test -entry PHOENIX --kraken2db $PATH_TO_DB
    ```

Note that this command clones (downloading) the repo to `~/.nextflow/assets/cdcgov/phoenix`. See [wiki](https://github.com/CDCgov/phoenix/wiki/Dependencies-and-Install#run-phoenix) for how to clone and have the software downloaded to a different location. 

    > * The pipeline comes with config profiles called `docker` and `singularity` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
 

7. Start running your own analysis with a [samplesheet](https://github.com/cdcent/phoenix/wiki/Running-PHoeNIx#samplesheet-input)!

    ```console
    nextflow run cdcgov/phoenix -r v1.0.0 -profile <singularity/docker/custom> -entry PHOENIX --input <path_to_samplesheet.csv> --kraken2db $PATH_TO_DB
    ```
    
# CDCgov GitHub Organization Open Source Project

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Access Request, Repo Creation Request

* [CDC GitHub Open Project Request Form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) _[Requires a CDC Office365 login, if you do not have a CDC Office365 please ask a friend who does to submit the request on your behalf. If you're looking for access to the CDCEnt private organization, please use the [GitHub Enterprise Cloud Access Request form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUQjVJVDlKS1c0SlhQSUxLNVBaOEZCNUczVS4u).]_

## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

## Overview

Describe the purpose of your project. Add additional sections as necessary to help collaborators and potential collaborators understand and use your project.
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
