Development
===========

Gitlab Continuous Integration (CI)
----------------------------------

Gitlab's feature to support basic continuous integration (CI_) is used to perform integration test on each push to
the gitlab repository. The setup is defined in the file `.gitlab-ci.yml`

To setup a new test runner, run the script `tests/gitlab_CI_docker/build_enpt_testsuite_image.sh` on the box which
is supposed to run the test. Tests are run in docker images which are build by this script. Root access to the machine
is therefore needed. The setup of the test docker image is defined in the file:
`tests/gitlab_CI_docker/context/enpt_ci.docker`


.. _CI: https://about.gitlab.com/features/gitlab-ci-cd/
