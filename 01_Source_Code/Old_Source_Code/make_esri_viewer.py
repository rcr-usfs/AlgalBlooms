"""
   Copyright 2021 Ian Housman, RedCastle Resources Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
#Script make a basic ESRI viewer using an existing template and geojson
####################################################################################################
import json,os
##################################################
#Function to set up a viewer using ArcMap API for JavaScript
#Requires a template file
#The replace_dict is intended to replace any key with the value
# E.g. replace_dict = {'{REPLACE_THIS}','hello world'} would replace anywhere in the teplate with {REPLACE_THIS} with hello world
def setup_viewer(output_viewer_name, replace_dict, viewer_template):
  o = open(viewer_template,'r')
  viewer_lines = o.read()
  o.close()

  for k in replace_dict.keys():
    viewer_lines = viewer_lines.replace(k,replace_dict[k])

  o = open(output_viewer_name,'w')
  o.write(viewer_lines)
  o.close()
