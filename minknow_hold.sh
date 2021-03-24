#!/bin/bash

# Mark all MinKNOW related packages as auto
sudo apt-mark auto minknow-core-minion-nc
sudo apt-mark auto ont-bream4-minion
sudo apt-mark auto ont-configuration-customer-minion
sudo apt-mark auto ont-kingfisher-ui-minion
sudo apt-mark auto ont-vbz-hdf-plugin
sudo apt-mark auto minknow-nc


#!/bin/bash

# Mark all MinKNOW related packages as hold
sudo apt-mark hold minknow-core-minion-nc
sudo apt-mark hold ont-bream4-minion
sudo apt-mark hold ont-configuration-customer-minion
sudo apt-mark hold ont-kingfisher-ui-minion
sudo apt-mark hold ont-vbz-hdf-plugin
sudo apt-mark hold minknow-nc

