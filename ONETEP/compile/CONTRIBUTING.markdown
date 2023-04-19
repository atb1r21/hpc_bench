# CONTRIBUTING
## Guidelines for developing for ONETEP

James C. Womack

Nicholas D. M. Hine

Version 1.3 (12/09/18)

These guidelines are adapted from guidelines for contributing to the CASTEP project, kindly shared by the CASTEP developers.

### Summary
The ONETEP project is developed within a Git repository hosted on Bitbucket.

The project is hosted as a private repository, owned by the "ONETEP" team. The team has an academic account, which allows unlimited contributors and provides some free build minutes on Bitbucket's "Pipelines" service.

We have decided to adopt a workflow based on forking of an official ONETEP repository. To contribute to the ONETEP project:

1. Create a private fork of the [official ONETEP repository](https://bitbucket.org/onetep/onetep).
2. Clone this repository to your local machine.
3. Make changes to the repository.
4. When you are ready to contribute your changes, create a [pull request](https://www.atlassian.com/git/tutorials/making-a-pull-request).

Your proposed changes will be reviewed by other developers and, if accepted, merged into the official repository.

### Signing up for Bitbucket

* To contribute to Bitbucket, you need a Bitbucket account.
* If you already have a Bitbucket account, we recommend that you [request an academic license](https://www.atlassian.com/software/views/bitbucket-academic-license) (if you are eligible).
* If you do not already have a Bitbucket account, we recommend you [sign up](https://id.atlassian.com/signup) using an academic e-mail address, as this should automatically put your account on an academic plan.
* In order to connect to Bitbucket using SSH (recommended), you will need to [add an SSH key to your Bitbucket account](https://confluence.atlassian.com/bitbucket/set-up-an-ssh-key-728138079.html).

### Gaining access to the ONETEP repository on Bitbucket
To access the ONETEP repository, you need to become a member of the [ONETEP team](https://bitbucket.org/onetep). Please contact one of the team administrators with your Bitbucket username to request this.

Current administrators include:

* Jacek Dziedzic
* Nick Hine
* Chris Skylaris
* Joseph Prentice

You will be added to the "Contributors" (onetep:contributors) group. This gives you read access for the official repository and allows you to make private forks of this repository. Members of the contributors group cannot directly push changes to the official repository, but must contribute these via pull requests.

### Forking
To create a private fork of the [official ONETEP repository](https://bitbucket.org/onetep/onetep):

1. Log in to [Bitbucket](https://bitbucket.org).
2. Open the  [official ONETEP repository](https://bitbucket.org/onetep/onetep) project page in a web browser.
3. Open the "+" menu from the sidebar on the left of the page and select "Fork this repository".
4. In the page that opens, set the name of your fork, ideally by appending your initials, e.g. `onetep_JCW`, and then select "Fork repository".
5. You now have your own private fork of the official ONETEP repository.
6. You can now clone your fork to your local machine, commit changes and push these to the repository on Bitbucket without directly affecting the official repository.
7. **Optional** To allow other ONETEP developers to view your private fork, we recommend that you give the "onetep:contributors" group read access. To do this:
    * Go to the Bitbucket project page for your private fork in a web browser.
    * Select "Settings" from the sidebar menu, then "User and group access".
    * Use the drop down menu under "Groups" to select "Contributors (onetep:contributors)" and add read access.

### Development within a fork

* Once you have created your private fork, you can clone this to your local machine and start making changes, e.g. `git clone git@bitbucket.org:<username>/<repo_name>.git`, where `<username>` is your bitbucket username and `<repo_name>` is the name you gave your private fork.
* Remember to commit your changes regularly and push these to your fork on Bitbucket so that they are backed up, e.g. `git push origin <branch_name>` where `<branch_name>` is the name of the git branch you are working on.
* You can choose how you develop in your private fork as the branches, tags etc. that you create in the fork do not affect the official repository.
* Keep your fork up-to-date with changes in the official repository by regularly merging with the official repository, i.e.
    * Add the official repository as a remote to your local repository: `git remote add bitbucket_official git@bitbucket.org:onetep/onetep.git`
    * Fetch changes from the official repository: `git fetch bitbucket_official`
    * Merge changes from the official repository into your currently checked-out branch: `git merge bitbucket_official/master`
* You can apply bugfixes to release branches in your private fork, too, i.e.
    * Fetch changes from the official repository: `git fetch bitbucket_official`
    * Check out a release branch: `git checkout academic_release_v5.0`
    * Commit your bugfix to the branch and push to your changes back to Bitbucket.

The above text describes a basic Git workflow within a private fork. If you are not familiar with Git, or source code version control in general, it may be worthwhile to spend some time working through a tutorial or guide before proceeding with any serious ONETEP development. There are numerous resources available online. Here is a small selection:

* [Atlassian Bitbucket Git tutorials](https://www.atlassian.com/git/tutorials)
* [Official Git (git-scm.com) documentation](https://git-scm.com/doc)
* [Software carpentry: "Version control with Git"](https://swcarpentry.github.io/git-novice/)

More experienced Git users may prefer to use a different process. You may use whatever workflow you like in your own private fork, since your changes will be selectively applied to the official repository via the pull request process.


### Creating a pull request

* When your changes are ready contribute to the official ONETEP repository, you need to create a pull request via the Bitbucket web interface.
* **Before you make your pull request**:
    * Make sure that your changes satisfy the code quality requirements (see below), or your pull request will be rejected.
    * Make sure that your code is up-to-date with official ONETEP repository -- use the instructions above to merge changes into the branch you are working on in your fork.
* Once you are sure your code satisfies the requirements and is synchronized with the official repository, open the Bitbucket project page for your private fork in a web browser.
* Select "Pull requests" from the side menu, then "Create pull request".
* You will now be presented with an interface which allows you to configure the pull request:
    * On the left, select the branch on your private fork that has new commits you would like to merge into the official repository.
    * On the right, select the branch on the official repository where these commits should be applied.
    * For feature development, you should be merging into the `master` branch of the official repository.
    * For bugfixes to releases, you should be merging into the corresponding release branch of the official repository, e.g. `academic_release_v5.0`.
* Enter a descriptive title for your pull request in the "Title" field.
* Enter a description of the changes being contributed in the "Description" field.
* Select appropriate reviewers (from the "Developers" group in the ONETEP team) for your pull request in the "Reviewers" field.
* Select "Create pull request".

Now your pull request will be reviewed by members of the "Developers" group in the ONETEP team. If the pull request is accepted by the reviewers, it will be merged into the official repository.

If there are issues with your pull request, you may be asked to make changes. There is no need to create a new pull request. You can apply the requested changes directly to the branch in your private fork for which the pull request was issued and they will be automatically added your pull request.

### Code quality requirements

#### Pre-pull-request checks
The ONETEP source code is distributed with a number of useful scripts for checking the correctness of source code. Before creating a pull request, you should ensure that these tests all pass for the branch you want to merge in your private fork:

1. **Check module dependencies**: Run `make checkdeps` in the root directory of the repository.
2. **Check input keywords**: Run `make checkesdf` in the root directory of the repository.
3. **Check whitespace**: Run `./utils/check_whitespace` in the root directory of the repository. This script can automatically fix whitespace issues -- see the output of `./utils/check_whitespace -h` for details.

You should also run the full QC test suite before creating a pull request. This can be done by running `make ARCH=<ONETEP_ARCH>` in the `tests` directory of the repository (`<ONETEP_ARCH>` is the name of the configuration file you used to compile ONETEP, e.g. `RH7.intel17.omp.scalapack`). Key points about running the QC test suite before creating a pull request:

* You should run the test suite with ONETEP compiled with a recent version of a widely used Fortran compiler, e.g. GFortran, Intel Fortran.
* The tests should be run in parallel (ONETEP must be compiled with MPI support), ideally also making use of OpenMP threading (set the environment variable `OMP_NUM_THREADS` to a value >1).
* No tests should fail, though warnings are acceptable.

If tests fail, please investigate this further before creating your pull request:

* Check out the master branch from the official repository (ensuring that it is fully updated by pulling from it) and run the failing tests on this.
* If the relevant tests do not fail upstream (i.e. for the master branch on the official repository), you should fix the code in your fork before creating your pull request.
* If the relevant tests also fail upstream, please check the [list of issues in the official repository](https://bitbucket.org/onetep/onetep/issues): there should be an open issue with relevant information. If such an issue has not been reported, please create a new issue explaining which tests are failing and any details about how the failures occur.

If the tests are not failing upstream, the problem is likely to be part of the changes you have implemented and you will be the best person to find it and fix it. If you need help, and want to discuss the issue with other developers, you may create an issue in your fork (not in the official repository) and draw their attention to it via email or on the [ONETEP development Slack](https://onetep-devel.slack.com/).

If your tests fail in a way in which the official master branch is known to fail, your pull request may be acceptable -- you can go ahead creating the pull request, but make sure that the text of the pull request includes a list of the tests that fail and why you think that this is not a problem.

Finally, you should ensure that the ONETEP version number in your `src/onetep.F90` file is updated appropriately (see below for a description of ONETEP's version numbering system). If multiple pull requests are open, or it takes time to merge your pull request, it may be necessary to update the version number again during the pull request code review. This is to ensure that the version number increments linearly in the official repository's `master` branch.

#### Style conventions
In order to keep the ONETEP codebase relatively consistent and readable, we have adopted the following conventions:

* Each line of code should be 80 columns or less in length, as this makes comparison of two versions of code side-by-side easier (some older code does not adhere to this, but all new code should).
* If a line exceeds 80 columns in length, use continuation characters (`&`) to break it into chunks of ≤80 characters (remember to inclide an extra '&' at the beginning of the following line when breaking in the middle of a string literal).
* Always uses spaces for indentation, never tabs.
* The blocks in the `do` loop, `select case` (also `select type`) and `if` constructs should be indented by 3 spaces.
* The contents of subroutines and functions should be indented by 2 spaces.
* Line continuations should be indented by 5 spaces for continued lines.
* Ensure there is no trailing whitespace before you commit (you can use the `./utils/check_whitespace` script described above to check this).
* In general, adding or removing blank lines should be avoided in core modules, as these changes will appear in the commit history.

#### New functionality
If you add new procedures or significantly change existing procedures, **you must create or update the documentation in the source code**. Examples of how to  document procedures can be found in the template module `template_mod.F90` in the `doc` directory. The key components of the documentation of procedures are:

* A human-readable description of the new functionality.
* A list of the arguments and a description of their meaning.
* The author(s) and a changelog describing significant modifications.

When adding new functionality which does not fit into other modules, it may be necessary to create a new source file containing a new module. Note that procedures and variables should always be encapsulated in modules, not 'bare' in a source file.

Before creating a new module, you should consider carefully whether your new functionality fits within the framework of an existing module, or is generic enough to be part of a multi-purpose module, such as `utils` or `services`. If a new module is needed to encapsulate some new functionality, then you should follow the following guidelines:

* Give your module a name which indicates the functionality it contains. If unsure, consult a more experienced developer to discuss an appropriate name.
* The filename for the module should have the form `<module_name>_mod.F90` where `<module_name>` is the name you have given the module.
* By default, variables and procedures in your module should be private (i.e. they should have the `private` attribute).
* Global variables (private or public variables declared at the level of the module itself rather than within its routines) constitute "hidden state", which tend to make the behaviour of a routine undesirably dependent on more than just the arguments it is called with. Sometimes these are unavoidable, and there are instances of them in the code. However, they should be minimised as much as possible. Think carefully before declaring any module-level global variables. More experienced developers may be able to suggest ways to encapsulate data inside arguments to routines such that they do not constitute "hidden state".
* Variables and procedures which do have to be public (accessible outside the module) should be explicitly specified (i.e. they should have the `public` attribute).
* In general, public variable and procedure names should be prepended by a standard prefix (typically the module name, or a shortened version of the name).

It is recommended that you make a copy of `./doc/template_mod.F90` and use this as a starting point for your new module, as this will make following the above guidelines easier.

#### Version numbering
There are four parts to the version number of a development version, and three parts to that of a release version. The first version number is only very rarely incremented by a collective decision of the main authors of the code (ODG). New major versions are released around every 6-12 months and are indicated by incrementing the second number in the full version number (e.g. "2" in "4.2").

The major version number (second number in full version number) indicates whether the associated source code is a release version or a development version:

* Release versions (which are distributed to users) have an _even_ second number.
* Development versions (which are under active development) have an _odd_ second number.

Within a series with the same second version number, successive versions (indicated by the third number in the full version number) should be compilable and complete with respect to a given new feature.

Release versions have three numbers (e.g. 4.2.3), while development versions have four numbers (e.g. 4.2.3.4).

The version number of a development version should be incremented prior to creating a pull request to be merged into the official repository's master branch. For minor changes (e.g. a bugfix or minor change to existing code), the fourth number should be incremented, while major changes (e.g. a new module or overhaul of existing functionality) should increment the third number (with the fourth number being reset to zero).

Bugfixes to a release version (merged into the corresponding release branch on the official repository, e.g. `academic_release_v5.0`) should increment the last (third) number in the full release version number.

**At any given time, there is a development version and a release version differing in their second version number by 1.**

    4.0.0   <-- first release of v4
    4.0.1   <-- bugfix to v4 in git branch for release

    4.1.0.0 <-- first development version of v4 (initially same as 4.0.0)
    4.1.0.1 <-- minor development work
    4.1.1.0 <-- significant development work
    4.2 RC3 <-- release candidate 3 for v4.2
    4.2.0   <-- next release version
    4.3.0.0 <-- next development version (initially same as 4.2.0).

    5.0.0   <-- first release of v5

### Notes for members of the "Developers" group
#### Preventing accidental pushes to the official repository

* Members of the "Developers" group (onetep:developers) in the ONETEP team, have write access to the official repository.
* Developers may want to take steps avoid accidentally pushing work to the official repository if they have added this as a remote to their private fork.
* This can be achieved by setting the push address for the remote to an unresolvable URL, i.e. `git remote set-url --push bitbucket_official DISABLE`.

### History

The earliest versions of the code, dating back to before 2005, were committed to a revision control system based on CVS. These files are still available on the TCM filesystem at /u/fs1/onestep/CVS_REPOSITORY. In around 2010 we moved to Subversion, using a repository still hosted on the TCM filesystem in Cambridge. The SVN history was migrated to Bitbucket and can still be browsed within the current Git repository (though attribution to authors is often not correctly recorded).

In June 2018, the ONETEP project was migrated from a Subversion repository to a Git repository. The hosting of source code was simultaneously moved from cvs.tcm.phy.cam.ac.uk to Bitbucket.

