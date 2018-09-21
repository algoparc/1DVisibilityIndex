#!/bin/bash

#Copyright 2016-2018 Ben Karsin, Nodari Sitchinava
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


echo "AoT sequential runtime"
for i in {16..20}
do
echo -n "Size 2^$i: "
../src/visAoT $i 1
done

echo "AoT 4 threads"
for i in {16..20}
do
echo -n "Size 2^$i: "
../src/visAoT $i 4
done

echo "AoT 8 threads"
for i in {16..20}
do
echo -n "Size 2^$i: "
../src/visAoT $i 8
done
