#!/usr/bin/env python
#
# Copyright 2009 Facebook
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may obtain
# a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.



import pymysql
pymysql.install_as_MySQLdb()
import pymysql as MySQLdb
import mysql.connector
import time
from .utils import get_config

class DB(object):

    def __init__(self):
        self._last_use_time = time.time()
        self.max_idle_time = float(1800)
        self._db = None
        args = {}
        config = get_config()
        args["user"] = config["database_user"]
        args["passwd"] = config["database_password"]
        args["db"] = config["database_name"]
        args["host"] = config["database_server"]
        args["port"] =  3306
        args["unix_socket"] = config["unix_socket"]
        self.host = "{0}:{1}".format(args['host'], args['port'])

        self._db_args = args


        """try:
            self.reconnect()
        except Exception:
            raise"""

    def __del__(self):
        self.close()

    def close(self):
        """Closes this database connection."""
        if getattr(self, "_db", None) is not None:
            self._db.close()
            self._db = None

    def reconnect(self):
        """Closes the existing database connection and re-opens it."""
        self.close()
        self._db = MySQLdb.connect(**self._db_args)

    def iter(self, query, params):
        """Returns an iterator for the given query and parameters."""
        self._ensure_connected()
        cursor = MySQLdb.cursors.SSCursor(self._db)
        try:
            self._execute(cursor, query, params)
            column_names = [d[0] for d in cursor.description]
            for row in cursor:
                yield row
        finally:
            cursor.close()

    def query(self, query, params=None):
        '''
        cursor = self._cursor()
        try:
            self._execute(cursor, query, params)
            result = cursor.fetchall()
            return result.
        except:
            self.reconnect()
            self.query(query, params)
        finally:
            cursor.close()
        '''
        #Potential fix for connection error from github copilot
        cursor = self._cursor()
        try:
            result = self._execute(query, params)
            return result
        except MySQLdb.OperationalError as e:
            print(("OperationalError: %s" % e))
            self.reconnect()
            return self.query(query, params)
        finally:
            if cursor:
                try:
                    cursor.close()
                except MySQLdb.OperationalError as e:
                    print(("Error closing cursor: %s" % e))


    def get(self, query, params=None):
        """Returns the first row returned for the given query."""
        rows = self._execute(query, None)
        if not rows:
            return None
        elif len(rows) > 1:
            raise Exception("Multiple rows returned for Database.get() query")
        else:
            return rows[0]

    def _ensure_connected(self):
        # Mysql by default closes client connections that are idle for
        # 8 hours, but the client library does not report this fact until
        # you try to perform a query and it fails.  Protect against this
        # case by preemptively closing and reopening the connection
        # if it has been idle for too long (7 hours by default).
        if (self._db is None or (time.time() - self._last_use_time > self.max_idle_time)):
            self.reconnect()
        self._last_use_time = time.time()

    def _cursor(self):
        self._ensure_connected()
        return self._db.cursor()

    def _execute(self, query, params):
        try:
            with mysql.connector.connect(**self._db_args) as connection:
                with connection.cursor() as cursor:
                    cursor.execute(query)
                    results =  cursor.fetchall()

        except MySQLdb.OperationalError:
            print(("Error connecting to MySQL on %s" % self.host))
            self.close()
            raise StarteratorError("Error connecting to database! Please enter correct login credentials in Preferences menu.")

        return results

class Row(dict):
    """A dict that allows for object-like property access syntax."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)


def get_db():
    database = DB()
    print("here")
    return database