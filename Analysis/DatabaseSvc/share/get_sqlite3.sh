#!/bin/sh

SQLITE_URL="http://www.sqlite.org/sqlite-amalgamation-3070500.zip"
LOCAL_ZIP="sqlite3.zip"
if [ ! -f ${LOCAL_ZIP} ]; then
  wget ${SQLITE_URL} -O ${LOCAL_ZIP}
fi
unzip -u -j ${LOCAL_ZIP} "sqlite-amalgamation-*/sqlite3.*"

