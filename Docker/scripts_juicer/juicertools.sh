#!/bin/bash

command="java -Xms512m -Xmx64384m -jar /opt/juicer/scripts/common/juicer_tools.jar $@"
echo $command
eval $command
