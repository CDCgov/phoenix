# CDC GitHub Practices for Open Source Projects

**The [CDCGov organization on GitHub](https://github.com/CDCgov) is designated for use by CDC programs to publish open source code.** This is a set of practices to help programs release secure and compliant open source projects successfully. If you are interested in using GitHub for non-open source projects, please see information on our [enterprise organization](#cdc-enterprise).

We designed these practices to be straightforward and helpful, and we [accept feedback](#support-and-feedback) from the community on updating them. For [Required Practices](#required-practices), Projects that don't adhere to the [Required Practices](#required-practices) could be subject to [archival or removal](#non-compliance-procedure).

## Getting Started

Before you can publish your project, you must request access to be added to the CDCgov organization. Complete these steps:

1. Review the [Rules of Behavior](rules_of_behavior.md).
2. Confirm your [Github profile is setup](#profile-setup) properly.
3. Complete the [project request form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u).
   * This will require your CDC login, so if you don't have a login, ask someone to request on your behalf, or [get in touch](#support-and-feedback).

You should receive an email or notification when you are given access and your first repository should be setup for you. For subsequent projects, you will be able to create a repository in the organization using Github's interface. The [template repository](https://github.com/CDCgov/template) is maintained and an easy way to quick start your repository that complies with the guidelines. Once this is completed you're ready to follow the required guidelines to publish code.

## Required Practices

You must follow these practices before you publish real code into your repository.

* [ ] **Get Clearance.** Always obtain clearance from your organization prior to setting up and publishing a repository.
  * GitHub is a third party service used by CDC to collaborate with the public. Official CDC health messages will always be distributed through www.cdc.gov and through appropriate channels, so make sure to plan your project along with your official public health program on cdc.gov.
* [ ] **Naming.** Set a meaningful project name and short description for your project. The form to do this is in your repositories settings.
  * [ ] Add [topics](https://help.github.com/en/github/administering-a-repository/classifying-your-repository-with-topics) to improve discovery and use of your project. For AI-related projects, the [Code.gov Implementation Guidance to Federal Agencies Regarding Enterprise Data and Source Code Inventories](https://code.gov/federal-agencies/compliance/inventory-code) must be followed when setting topics.
* [ ] **Create a README.** Add a `README.md` file at the root with the following:
  * An overview of your project, including the purpose, goals and the team responsible.
  * A description of your development process in the `README.md` file. If your project is no longer active, mark it as [archived](https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/archiving-repositories).
  * Include the following notice sections. You can modify the verbiage and adapt as necessary based on your program need.
    * [ ] [Public Domain Standard Notice](https://github.com/CDCgov/template#public-domain-standard-notice)
    * [ ] [License Standard Notice](https://github.com/CDCgov/template#license-standard-notice)
    * [ ] [Privacy Standard Notice](https://github.com/CDCgov/template#privacy-standard-notice)
    * [ ] [Contributing Standard Notice](https://github.com/CDCgov/template#contributing-standard-notice)
    * [ ] [Records Management Standard Notice](https://github.com/CDCgov/template#records-management-standard-notice)
    * [ ] [Additional Standard Notices](https://github.com/CDCgov/template#additional-standard-notices)
* [ ] **Choose a license.** Assign an open source license based on program need.
  * If you need help choosing a license, please review [this article](https://www.philab.cdc.gov/index.php/2012/03/27/open-source-development-for-public-health-informatics/), refer to existing CDCgov projects, or ask for consultation support in choosing a license.
* [ ] **Security scanning and review.**
  * **This is the final step before publishing and the most critical.**
  * All source code used within CDC systems must comply with all cybersecurity processes prior to production use, including static and dynamic scanning. The same applies to code published as open source.
    * If you are unsure about compliance, reach out to your organization's security officers.
  * Never commit sensitive information, including usernames, passwords, tokens, PII, PHI. To automate this, you can integrate pre-commit tools like [Clouseau](https://github.com/cfpb/clouseau) to systematically review material before committing.
    * Make sure that the commit history of your Github repository also doesn't have these things. In many cases it's easier to start a new repository and push up the code that has all sensitive information removed as the first commit.
  * Enable [GitHub automated security alerts](https://help.github.com/en/github/managing-security-vulnerabilities/about-security-alerts-for-vulnerable-dependencies) and configure notification for the repo admin to see.
* [ ] **Setup your profile.** [Active project committers need to add profile info to help collaboration.](#profile-setup)
  * [ ] **Two-factor authentication (2FA).** [Project admins must secure their account with two-factor-authentication.](https://docs.github.com/en/enterprise-server@2.21/github/authenticating-to-github/securing-your-account-with-two-factor-authentication-2fa)
* [ ] **Maintain your repository.** Once your repository is published, you must do the following to remain in compliance:
  * [ ] **Respond to critical security issues and communication from administrators.** Ignoring security issues or not responding to communication from administrators can result in [archiving or removal](#non-compliance-procedure).
  * [ ] **Archive old projects.** If you're no longer updating the project or have moved it's location, update your `README.md` file to let users know and [archive the repository](https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/archiving-repositories).

## Recommended Practices

Optional improvements to make your open source project more successful.

* [ ] Establish pull request templates to make it easier for contributors to send pull requests. For example [SDP-V has a checklist for each PR to match their development practices.](https://github.com/CDCgov/SDP-Vocabulary-Service/blob/master/.github/PULL_REQUEST_TEMPLATE)
* [ ] Agree on project conventions and include them in your `README.md` file. Depending on what type of project, this includes folder structure for data, linters, editor configuration (eg, [MicrobeTrace's .editorconfig](https://github.com/CDCgov/MicrobeTrace/blob/master/.editorconfig)). This will help improve the quality of your project and make it easier for others to contribute to your project.
* [ ] Add support and community procedures. CDC does not provide warranty or official support for open source projects, but describing how you would like questions and issues will assist users of your project. If you use a wiki, or project board, or package manager, describe and link to that. Official contribution steps will make it easier for people outside of CDC to contribute to your project.
* [ ] Include references to publications, presentations, and sites featuring your project.
* [ ] Add an entry to [open.cdc.gov](https://open.cdc.gov) to the [data](https://open.cdc.gov/data.html), [code](https://open.cdc.gov/code.html), [api](https://open.cdc.gov/apis.html), or [event](https://open.cdc.gov/events.html) page to help people find your project on cdc.gov
* [ ] Add versions and tags describing major releases and milestones. For example, [open.cdc.gov's releases each time a new version is published to the web site](https://github.com/CDCgov/opencdc/releases/tag/v1.0.9) or [geneflow's changelog](https://github.com/CDCgov/geneflow/blob/master/CHANGELOG.md).
* [ ] Follow [Semantic Versioning 2.0.0](https://semver.org/) when creating versions for your project.
* [ ] Describe and test reproducible practices to install and build your project. For example, [injury_autocoding's code section on running the project's scripts](https://github.com/cdcai/injury_autocoding#code)).
* [ ] Recognize contributors and existing resources that have helped the project. For example, [fdns-ms-hl7-utils' AUTHORS file](https://github.com/CDCgov/fdns-ms-hl7-utils/blob/master/AUTHORS).
* [ ] Automate build and test procedures to reduce the effort of outside contributors to send pull requests (eg, [Travis CI](https://travis-ci.org/), [Circle CI](https://circleci.com/), [GitHub Actions](https://help.github.com/en/actions))
* [ ] [Appropriately gather metrics](https://opensource.guide/metrics/) on how your project is used and incorporate this into your feature planning process.
* [ ] [Incorporate documentation into your development cycle](https://github.com/GSA/code-gov-open-source-toolkit/blob/master/toolkit_docs/documentation.md), and where possible, automate the generation of documentation so it is more likely to be up to date and useful to people interested in your project.

## Guidance

### Support and Feedback

If you need additional support with your setting up project, or have any feedback or ideas about this guidance please [open an issue](https://github.com/CDCgov/template/issues) or send an email to [data@cdc.gov](mailto:data@cdc.gov). We also accept pull requests if you want to directly edit the guidance.

### Non-Compliance Procedure

Projects in this organization are reviewed occasionally for compliance with the [Required Practices](#required-practices). If your project is found to not be in compliance, you will be contacted by administrators to help bring your project into compliance. Projects that do not respond or that habitually fail to meet these practices will be archived or removed from the organization, depending on severity.

### Profile Setup

Please make sure your profile is set up properly to help us work better together. Specifically, keep your profile up to date with:

* **Name:** Your first and last name.
* **Company:** Your government agency or contracting company. (If you also use GitHub for personal projects, consider specifying “CDC (work) + personal projects” to make it clear that some of your GitHub projects may be personal in nature.)
* **Location:** Your primary work location (city, state).
* **Photo:** A headshot photo, or an appropriate image that is unique to you.

If you admin any projects, make sure to [secure your account with two-factor authentication (2FA)](https://docs.github.com/en/enterprise-server@2.21/github/authenticating-to-github/securing-your-account-with-two-factor-authentication-2fa). Although you probably already did this because you are smart.

### Open Source Checklist

So you've decided to set up an open source project at CDC. Here are the steps to do that, in the most common order.

* [ ] Create a new project using the [template repo](https://github.com/CDCgov/template).
* [ ] Update your readme.md following the [CDC GitHub Practices for Open Source Projects](https://github.com/CDCgov/template/blob/master/open_practices.md)
* [ ] Choose a license. Most projects are ASL2, but license should meet public health program need. See <https://www.philab.cdc.gov/index.php/2012/03/27/open-source-development-for-public-health-informatics/> for more info on choosing a license.
* [ ] Remove all sensitive info.
* [ ] Talk with your ADI, ADS, and ISSO for review and clearance.
* [ ] After approval, create a GitHub user.
* [ ] Fill out the [Request a Repo form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) for a new repo on [CDCGov](https://github.com/cdcgov) or [CDCai](https://github.com/cdcai).
* [ ] When you get an email or push alert that your repo is ready, push to GitHub
* [ ] Add an entry in [open.cdc.gov](https://open.cdc.gov) on their [code page](https://open.cdc.gov/code.html) to officially be linked from cdc.gov. This helps users find and use your project.
* [ ] Keep your project up to date, when you're finished flag it as [archived](https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/archiving-repositories).

_This checklist was adapted from the CDC IT Guard Rail and put here to help people who don't have access to the intranet._

### CDC Enterprise

Our [CDCent](https://github.com/cdcent/) organization is used for private, non-public projects so only CDC staff and approved outside collaborators work on these projects, you can request access through the [GitHub Enterprise Cloud form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUQjVJVDlKS1c0SlhQSUxLNVBaOEZCNUczVS4u).

### Reference Links

These are helpful links from across the Federal Government regarding open sourcing code.

* [CFPB Open Tech](https://cfpb.github.io/)
* [TTS Engineering Practices Guide](https://engineering.18f.gov/)
* [18F Open Source Policy](https://github.com/18F/open-source-policy) and [Practicing our open source policy](https://github.com/18F/open-source-policy/blob/master/practice.md)
* [GitHub and Government: How agencies build software](https://government.github.com/)
* [code.gov](https://code.gov)
* [Federal Source Code and Open Source Toolkit](https://github.com/GSA/code-gov-open-source-toolkit)
* [Federal Source Code Policy (M-16-21)](https://sourcecode.cio.gov/)
* [openCDC](https://open.cdc.gov)
* [Digital Services Playbook](https://playbook.cio.gov/)
* [CDC/ATSDR Policy on Public Health Research and Nonresearch Data Management and Access](https://www.cdc.gov/maso/policy/policy385.pdf)
  * [CDC/ATSDR Policy on Releasing and Sharing Data](https://www.cdc.gov/maso/Policy/ReleasingData.pdf) (old version, but still a useful reference)
* [Clearance of Information Products Disseminated Outside CDC for Public Use](https://www.cdc.gov/os/policies/docs/CDC-GA-2005-06_Clearance_of_Information_Products_Disseminated_Outside_for_Public_Use.pdf)
* [Federal Source Code Toolkit](https://github.com/GSA/code-gov-open-source-toolkit)
