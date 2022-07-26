#!/bin/bash

command="java -Xms512m -Xmx32384m -jar /opt/juicer_tools.jar $@"
echo $command
eval $command
