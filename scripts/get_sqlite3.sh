#!/bin/sh

# SQLITE_URL="http://www.sqlite.org/sqlite-amalgamation-3070500.zip"
SQLITE_URL="https://www2.sqlite.org/2022/sqlite-amalgamation-3390200.zip"
LOCAL_ZIP="sqlite3-amalgamation.zip"
if [ ! -f ${LOCAL_ZIP} ]; then
  wget ${SQLITE_URL} -O ${LOCAL_ZIP}
fi
unzip -u -j ${LOCAL_ZIP} "sqlite-amalgamation-*/sqlite3.*"

