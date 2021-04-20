#!/bin/bash
#Mark all MinKNOW related packages as auto
sudo apt-mark unhold minion-nc
sudo apt-mark unhold minknow-core-minion-nc
sudo apt-mark unhold minknow-nc
sudo apt-mark unhold ont-bream4-minion
sudo apt-mark unhold ont-configuration-customer-minion
sudo apt-mark unhold ont-kingfisher-ui-minion
sudo apt-mark unhold ont-vbz-hdf-plugin

sudo apt-mark auto minion-nc
sudo apt-mark auto minknow-core-minion-nc
sudo apt-mark auto minknow-nc
sudo apt-mark auto ont-bream4-minion
sudo apt-mark auto ont-configuration-customer-minion
sudo apt-mark auto ont-kingfisher-ui-minion
sudo apt-mark auto ont-vbz-hdf-plugin

