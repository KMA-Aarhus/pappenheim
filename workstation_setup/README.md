# Workstation setup

### Install Ubuntu LTS


1. Download latest Ubuntu LTS (long term support) ISO image
2. Make a bootable Ubuntu USB stick.
3. Turn off the workstation.
4. Insert the bootable USB.
5. Turn on the workstation and press `F12` until boot menu opens up.
6. Choose the USB boot device from the boot menu to boot up with the Ubuntu installer.
7. In the Ubuntu installer menu, choose `Install Ubuntu`.
8. Language: English, keyboard layout: Danish.
7. In the **Updates and other software window** choose `minimal installation` and `Download updates while installing Ubuntu`.
   * Select **Install third-party software for graphics [...]**
9. In the **Installation type** window choose `Something else`.
10. In the **Partition window** erase the two disks
   * Delete all of the old partitions on SSD and HDD. This will delete the windows installation.
   * SSD
     * Choose add on free space. Choose EFI under Use as, with 600 mb in size.
     * Creat ROOT partition with remaining space (EXT4). Mount point: `/` ("Root")
   * HDD (or second larger SSD)
     * Create partition which maps to /home using all space (location of users home directories). Same as ROOT partition, but mount point is `/home`
   * Click `Install now`.
11. Set up a reasonable account. 
12. Print a sticker of computer name and put it on the bottom frame in the front.

NOTE: if black screen after restart, switch monitor cable ports.

### Update Ubuntu

1. Update Ubuntu packages by running the following command in a terminal:
   ```
   # Update the OS and packages
   sudo apt update
   sudo apt upgrade -y
   
   # Install latest versions useful software:
   sudo apt install git vim tree htop libreoffice filezilla tmux
   ```


### Install Nvidia and CUDA drivers




1. Check that the graphics driver is installed
  ```
  ubuntu-drivers devices
  ```

1B. If not, then install drivers and reboot:
  ```
  sudo ubuntu-drivers autoinstall
  ```
   
2. Make sure the CUDA drivers are installed as well

  ```
  sudo apt install nvcc
  ```

4. Check nvidia driver is working

  ```
  nvidia-smi
  ```

Notice: You may have to restart the computer between some of these steps. Who knows.



### Install MinION software and Guppy-gpu basecaller

1. Install "MinION Software" according to the instructions given on the minion download page
   https://community.nanoporetech.com/downloads
   
2. Now things become a bit fishy: The version of guppy you need is NOT the one that is stated on the download page.
   The way to install guppy is to find a compatible version of guppy in relation to the MinKNOW Software. See the table here: https://hackmd.io/@Miles/B1U-cOMyu#Guppy--MinKNOW-compatibility
   A nanopore forum user has been helpful giving the link to a patched version that works with cuda 11+ : https://community.nanoporetech.com/posts/minknow-and-guppy-version
   Hint: use `apt list minion-nc` to see the version of minion-nc
   
3. Then map the guppy-gpu installation to the MinKNOWN guppy folder:

   ```
   sudo /opt/ont/minknow/bin/config_editor --conf application --filename /opt/ont/minknow/conf/app_conf \
       --set guppy.server_executable="/path/to/guppy/bin/guppy_basecall_server" \
       --set guppy.client_executable="/path/to/guppy/bin/guppy_basecall_client" \
       --set guppy.gpu_calling=1 \
       --set guppy.num_threads=3 \
       --set guppy.ipc_threads=2
   ```

 Restart the minknow service or restart the computer, for minknow to register the changes. 
   
NOTE: If the installation worked correctly, `Guppy basecaller server` should be reserving a part of the GPU card and Guppy should be recognised by MinKnow. If not, MinKnow will not commence sequencing and give `Internal error` after the flow cell has reached sequencing temperature. This issue is usually caused by software incompatibility  and can often be solved by installing a different version of Guppy.

4. Install Miniconda3:
   
    ```
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    ```

5. Install rampart by following the [official instructions](https://github.com/artic-network/rampart/blob/master/docs/installation.md#install-from-conda).

    * After following the install instructions, and checking that the install is successfull, add a QnD alias to the shell-environment.
    But first, download this git repo which will be needed:
    ```
    mkdir -p ~/repos
    cd ~/repos
    git clone https://github.com/artic-network/artic-ncov2019.git
    ```
    
    Then, insert this into your ~/.bashrc
    ```
    alias start_rampart='{ $(sleep 3; firefox localhost:3000) & }; conda activate artic-rampart && rampart --protocol ~/repos/artic-ncov2019/rampart/ --clearAnnotated --basecalledPath'
    ```
    
6. Install the Pappenheim pipeline according to its [instructions](https://github.com/KMA-Aarhus/pappenheim#installation)
   * Add the following bash alias:
   ```
   echo 'source ~/pappenheim/scripts/bash_extensions.sh' >> ~/.bashrc && source ~/.bashrc
   
   
   ```
7. Run this command to get a link of the laboratory-instruction onto the desktop.
   ```
   ln -s ${HOME}/pappenheim/laboratory/SARS-CoV-2\ sekventeringsinstruks.pdf ~/Desktop
   
   
### Miscellaneous

Disable printer discovery: https://cirovladimir.wordpress.com/2019/02/11/ubuntu-18-04-disable-network-printer-auto-discovery/



Install miniconda 3

Install sublime text https://www.sublimetext.com/3

Install typora https://typora.io/





Consider holding the current version of minion/minknow. This will make sure that no incompatibilities can arise when apt does automatic updates in the background.
```
~/pappenheim/workstation_setup/minknow_hold.sh 

# Reverse above script with the following:
# ~/pappenheim/workstation_setup/minknow_auto.sh 



```
