#!/bin/bash

# Ensure you are on the main branch before building.
build_documentation=false

# This will force-push to override the existing 'documentation' branch and
# trigger a GitHub Actions workflow to update the website.
deploy_documentation=false

# Delete the 'documentation' branch in the local GitHub repository, delete the
# 'docs' folder, and return to the 'main' branch.
cleanup_git=false

if [[ $# -gt 3 ]]; then
  echo "Too many arguments."
  invalid_input=true
elif [[ $# != 0 ]]; then
  invalid_input=false
 
  if [[ $invalid_input == false ]]; then
    for param in "$@"; do
      if [[ $param == "--build" ]]; then
        if [[ $build_documentation == true ]]; then
          echo "Duplicate argument '--build'."
          invalid_input=true
        else
          build_documentation=true
        fi
      elif [[ $param == "--deploy" ]]; then
        if [[ $deploy_documentation == true ]]; then
          echo "Duplicate argument '--deploy'."
          invalid_input=true
        else
          deploy_documentation=true
        fi
      elif [[ $param == "--cleanup" ]]; then
        if [[ $cleanup_git == true ]]; then
          echo "Duplicate argument '--cleanup'."
          invalid_input=true
        else
          cleanup_git=true
        fi
      else
        echo "Unrecognized argument '${param}'."
        invalid_input=true
      fi
    done
  fi
else
  echo "No arguments found."
  invalid_input=true
fi

if [[ $invalid_input == true ]]; then
  echo "Usage: documentation.sh [--build] [--deploy] [--cleanup]"
  exit -1
fi

if [[ $build_documentation == true ]]; then
  git branch -D documentation
  git checkout -b documentation
  rm -rf docs && mkdir docs
  swift package --allow-writing-to-directory docs generate-documentation --target MM4 --transform-for-static-hosting --output-path docs --hosting-base-path MM4
  echo '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=https://philipturner.github.io/MM4/documentation/mm4">' > docs/index.html
fi

if [[ $deploy_documentation == true ]]; then
  git add "docs/*"
  git commit -m "Update docs."
  git push --force origin documentation
fi

if [[ $cleanup_git == true ]]; then
  git checkout main
  git branch -D documentation
  rm -rf docs
fi
