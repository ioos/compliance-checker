#!/usr/bin/env bash

set -euxo pipefail

RELEASE_BRANCH=master

bumpversion_build() {
  bump2version patch
}

bumpversion_release() {
  git fetch --prune origin "+refs/tags/*:refs/tags/*"
  bump2version patch
  VERSION=$(bump2version --list --tag --commit --allow-dirty release | grep -oP '^new_version=\K.*$')
  git push origin tag $VERSION
  git push origin HEAD:$RELEASE_BRANCH
}

main() {
  local mode=$1; shift

  if [ "x${mode}" == "xbuild" ]
  then
    bumpversion_build
  elif [ "x${mode}" == "xrelease" ]
  then
    bumpversion_release
  fi

  exit 0
}

main "$@"
