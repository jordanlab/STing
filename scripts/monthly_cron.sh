#!/bin/env bash

dir=`dirname "$0"`
home_path=`realpath "$dir"`
cd "$home_path"

date=`date +"%Y%m%d-%H%M"`


# TODO: only mlst typing schema supported for now

cat fetch_list |parallel ./db_util.py fetch -o $HOME/pubmlst_dbs/mlst -q {} -b 2>&1 > cron-${date}.log

