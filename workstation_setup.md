# Workstation setup

### Install Ubuntu 18

1. Download Ubuntu 18 ISO image
2. Make a bootable Ubuntu USB stick by following these instructions.
3. Turn off the workstation.
4. Insert the bootable USB.
5. Turn on the workstation and press `F12` until boot menu opens up.
6. Choose the USB boot device from the boot menu to boot up with the Ubuntu installer.
7. In the Ubuntu installer menu, choose `Install Ubuntu`.
8. Language: English, keyboard layout: Danish with Win keys
7. In the **Updates and other software window** choose `minimal installation` and `Download updates while installing Ubuntu`.
8. In the **Installation type** window choose `Something else`.
9. In the **Partition window** erase the two disks
   * Delete all of the old partitions on SSD and HDD. This will delete the windows installation.
   * SSD
     * Choose add on free space. Choose EFI under Use as, with 600 mb in size.
     * Creat ROOT partition with remaining space (EXT4). Mount point: `/` ("Root")
   * HDD (or second larger SSD)
     * Create partition which maps to /home using all space (location of users home directories). Same as ROOT partition, but mount point is `/home`
   * Click `Install now`.
10. Set up a reasonable account. 
11. Print a sticker of computer name and put it on the bottom frame in the front.

NOTE: if black screen after restart, switch monitor cable ports.

### Update Ubuntu

1. Update Ubuntu packages by running the following command in a terminal:
   ```
   sudo apt update
   sudo apt upgrade
   ```

### Install Nvidia and CUDA drivers

### Install Nvidia and CUDA drivers

1. Install nvcc and make by running following commands in terminal:

   ```
   nvcc
   ```


3. Check that the graphics driver is installed

   ```
   ubuntu-drivers devices
   ```

If not, then install drivers and reboot:

   ```
   sudo ubuntu-drivers autoinstall
   ```

4. Check nvidia driver is working

   ```
   nvidia-smi
   ```



**References**

- https://gist.github.com/bendangnuksung/981408031699e0ddc50a6f6fdcf185c2
- https://medium.com/@naomi.fridman/install-conda-tensorflow-gpu-and-keras-on-ubuntu-18-04-1b403e740e25

### Install MinION software and Guppy-gpu basecaller

1. Install "MinION Software" according to the instructions given on the minion download page
   https://community.nanoporetech.com/downloads
   
2. Now things become a bit fishy: The version of guppy you need is NOT the one that is stated on the download page.
   The way to install guppy is to find a compatible version of guppy in relation to the MinKNOW Software. See the table here: https://hackmd.io/@Miles/B1U-cOMyu#Guppy--MinKNOW-compatibility
   A nanopore forum user has been helpful giving the link to a patched version that works with cuda 11+ : https://community.nanoporetech.com/posts/minknow-and-guppy-version
   
3. Then map the guppy-gpu installation to the MinKNOWN guppy folder:

     ```
     sudo /opt/ont/minknow/bin/config_editor --conf application --filename /opt/ont/minknow/conf/app_conf \
         --set guppy.server_executable="/path/to/guppy/bin/guppy_basecall_server" \
         --set guppy.client_executable="/path/to/guppy/bin/guppy_basecaller" \
         --set guppy.gpu_calling=1 \
         --set guppy.num_threads=3 \
         --set guppy.ipc_threads=2
     ```

4. You may also modify the MinKNOW `user_conf` file to change default output folder

   * Open folder `/opt/ont/minknow/conf`

   * Open a terminal in the folder an run the following command:

     `sudo gedit user_conf`

   * Replace current content of `user_conf`with the following and save before exiting:

     ```
     {
         "user": {
             "cereal_class_version": 2,
             "log": {
                 "cereal_class_version": 0,
                 "trace_categories": ""
             },
             "output_dirs": {
                 "cereal_class_version": 2,
                 "base": {
                     "value0": "/"
                 },
                 "logs": {
                     "value0": "/var/log/minknow"
                 },
                 "protocol_output_pattern": "{protocol_group_id}/{sample_id}/{start_time}_{device_id}_{flow_cell_id}_{short_protocol_run_id}",
                 "intermediate": {
                     "value0": "/var/lib/minknow/data/intermediate"
                 },
                 "complete_intermediate_reads": {
                     "value0": "/var/lib/minknow/data/queued_reads"
                 },
                 "reads": {
                     "value0": "/home/auh-covid19/Desktop/covid19_analysis"
                 },
                 "reads_tmp": {
                     "value0": "/var/lib/minknow/data/reads/tmp"
                 },
                 "fallback_reads_path": {
                     "value0": "/var/lib/minknow/data/reads/fallback/"
                 },
                 "fastq_tmp": {
                     "value0": "/var/lib/minknow/data/reads/fastq_tmp"
                 },
                 "fallback_fastq_path": {
                     "value0": "/var/lib/minknow/data/reads/fastq_fallback/"
                 },
                 "ping_queue": {
                     "value0": "pings"
                 }
             },
             "proxy": {
                 "cereal_class_version": 0,
                 "use_system_settings": true,
                 "auto_detect": true,
                 "auto_config_script": "",
                 "https_proxy": "",
                 "proxy_bypass": ""
             }
         }
     }
     ```

   * Restart workstation for changes to take effect.
   
NOTE: If the installation worked correctly, `Guppy basecaller server` should be reserving a part of the GPU card and Guppy should be recognised by MinKnow. If not, MinKnow will not commence sequencing and give `Internal error` after the flow cell has reached sequencing temperature. This issue is usually caused by software incompatibility  and can often be solved by installing a different version of Guppy.
   
### Miscellaneous

Disable printer discovery: https://cirovladimir.wordpress.com/2019/02/11/ubuntu-18-04-disable-network-printer-auto-discovery/

Install useful software:
``` 
sudo apt install git vim
sudo apt install sublime-text # add channel first
sudo apt install typora # add add channel first

```
  