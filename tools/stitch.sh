#!/bin/bash

################################################################################
################################ BLUEBOTTLE-1.0 ################################
################################################################################
#
#   Copyright 2012 - 2014 Adam Sierakowski, The Johns Hopkins University
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
# 
#   Please contact the Johns Hopkins University to use Bluebottle for
#   commercial and/or for-profit applications.
################################################################################

######################
# stitch.sh: stitch together PNG ParaView images (called anim) into a
# cross-platform AVI called anim.avi.
#
# Copyright Adam Sierakowski, The Johns Hopkins University 2014
######################

if [[ -z $4 ]] ; then
# if not all four arguments are given, show configuration
  echo "stitch.sh usage: ./stitch.sh png_dir png_name avi_name fps"
else
# else read in the arguments
  DIR=$1/$2.%04d.png
  NAM=$2.%04d.png
  AVI=$1/../$3.avi
  FPS=$4
  echo "Confirm that you want to run stitch.sh with the following parameters:"
  echo
  echo "png_dir  = $DIR"
  echo "png_name = $NAM"
  echo "avi_name = $AVI"
  echo "fps      = $FPS"
  echo
  NRM=$(bc -l <<< "scale=2; 30/$FPS")
  PTS="setpts=$NRM*PTS"
  echo "Confirm (y/N): "
  read RUN
  if [[ -z $RUN ]] ; then
    echo "Aborting stitch.sh"
  else
    if [ $RUN == y ] || [ $RUN == Y ] ; then
      echo
      echo "Launching ffmpeg..."
      echo
      echo "================================================================================"
      echo
      ffmpeg -r 30 -qscale 4 -f image2 -i $DIR -vcodec msmpeg4v2 -r 30 -vf "$PTS" $AVI
      echo
      echo "================================================================================"
      echo
      echo "stitch.sh complete"
    else
      echo "Aborting stitch.sh"
    fi
  fi
fi

