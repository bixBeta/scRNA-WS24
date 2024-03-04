### Week 4. scanpy, docker, advanced topics

Most of this exercise will be done in a Jupyter notebook. The notebook will be run in a docker container. In this document, there are instructions for getting the docker container started, copying the jupyter notebook, and logging into your jupyter notebook. Once you have logged into jupyter, you can explore the jupyter notebooks, which will demonstrate single-cell analysis using scanpy.


#### 1. Connect to the server

Find your assigned server name on [this website](https://biohpc.cornell.edu/ww/machines.aspx?i=165). Instructions to connect to the server can be found on [this page](https://biohpc.cornell.edu/lab/doc/remote_access.pdf).

#### 2. Pull the docker image
```
docker1 pull biohpc/scrna2024
docker1 images
```
The first command downloads the image from dockerhub, if not already on the server. The second command lists all currently available images on the server, you should see `biohpc/scrna2024` on the list.

#### 3. Clone the workshop repository
Chose a working directory, 'cd' to it, and then download a copy of this repository from github. For working dir, you may chose a directory in `/workdir/$USER`, or in `/home2/$USER`. In this example I will use `/workdir/$USER/scrna_workshop`.
```
# first, make the directory if it does not already exist, then cd into it
mkdir -p /workdir/$USER/scrna_workshop
cd /workdir/$USER/scrna_workshop

# now use git clone to make a copy of workshop materials in this directory
git clone https://github.com/bixBeta/scRNA-WS24.git
```
<details>
  <summary>About choosing a working directory</summary>
  <p>We usually recommend doing your work in some subdirectory of <code>/workdir/$USER</code>. <code>/workdir/$USER</code> is local disk, so input/output operations will be fastest using this directory. However, you only have access to <code>/workdir/$USER</code> during your reservation. If you want to save anything, you will need to copy it to more permanent storage (such as your home directory) before the reservation ends.</p>
  <p> More recently, we also have the option of using <code>/home2/$USER</code> for computation. This is networked storage, but configured to be fast enough to support computations. It still may not be as fast as local disk. The benefit is that <code>/home2/$USER</code> is permanent storage, so you will not lose any data when your reservation ends. However, it is subject to a quota (usually 200Gb), unless you or your lab has purchased additional storage.</p>
</details>

### 4. Start a docker container
The container is a virtual machine, running the operating system stored in the `biohpc/scrna2024` image. The image is designed to start a jupyter notebook as soon as the container starts. The command will look like this (read on for some explanation):
```
docker1 run -p 8009:8888 -v `pwd`:/data \
    -v /workdir/sc_workshop_2024:/workdir/sc_workshop_2024:ro  \
    --rm -d biohpc/scrna2024 /data
```
<details>
<summary>Details about the <code>docker1 run</code> command options</summary>
  <ul>
<li> <tt>-p 8009:8888</tt>: maps the port 8009 on your assigned server to port 8888 inside the container. By default, jupyter uses port 8888, but this is just the port inside the virtual machine. To access it on your computer, you need to map this to a port on the host machine. By using this command, you should be able to access the jupyter server at http://cbsuXXXXX.biohpc.cornell.edu:8009. (replace cbsuXXXXX with your assigned machine). If you get an error about the port being busy, it may be being used by another workshop participant, so try another port between 8009-8039. Remember that 8016 is already being used for RStudio.</li>
<li> <code>-v `pwd`:/data</code>: This mounts your current directory (returned by the <code>`pwd`</code> command) to the directory <code>/data</code> inside the docker container. At BioHPC, the <docker1>docker1</docker1> version of docker restricts the directories you can mount; you should be able to mount any directory in <code>/workdir</code>, <code>/local</code>, or <code>/home2</code> owned by you. <code>docker1</code> also automatically mounts <code>/workdir/$USER</code> to <code>/workdir</code> inside the container by default.</li>
<li><code>-v /workdir/scrna_workshop_2024:/data/scrna_workshop_2024:ro</code>: This is an additional mount command (you can mount as many directories as you like); this time mounting <code>/workdir/scrna_workshop_2024</code> on the server to <code>/data/scrna_workshop_2024</code> inside the container. The <code>:ro</code> at the end means that it will be a read-only mount. <code>docker1</code> allows read-only mounting of directories owned by others, if the owner has granted access (using <code>docker1 access</code> command). By mounting this directory, you will be able to access the example data without making your own copy.</li>
<li> <code>--rm</code> just means to delete the container when it finishes running. It is not necessary, but helps keep the server clear of old zombie containers.</li>
<li> <code>-d</code> means to run the container in the background.</li>
<li> The <code>/data</code> at the end of the command is passed to the container, this directory will be used as jupyter's root directory.</li>
  </ul>
</details>

If successful, the command should ouput a long string of characters - this is the container ID. Usually only the first 12 characters of the ID are used. You can also verify the container is running with the command `docker1 ps`. This lists all running containers on the server.

If not successful, you may need to try a different port instead of 8009 (see detail section above for more info about choosing a port).

### 5. Log into your jupyter notebook
  First, you need to get the notebook URL, which includes a long authentication token. Copy your container ID from the previous section, and replace the XXXXXX in the following command with the container ID:
  ```
  docker1 exec XXXXXXXX jupyter lab list
  ```
This will return the address of your jupyter notebook. It will start with http://XXXXXX:8888. You need to change the XXXXX part to your server name (i.e., cbsumma01.biohpc.cornell.edu), and change 8888 to the port you have mapped to (i.e., 8009).  You can do this manually, or automate this with the following code:

```
# change port and containerID below as appropriate
port=8009 
containerID=XXXXXXX

# then run this command
docker1 exec $containerID jupyter lab list | tail -n 1 |
   awk -F '[ /]' -v host=$HOSTNAME -v port=$port '{print "http://"host":"port$4}'
```
The last command should output the URL of your jupyter notebook. You should be able to connect to it if you are on the Cornell campus network. If you are off campus, you will need to set up a tunnel (see instructions in details below).
<details>
  <summary>Connecting from outside the Cornell network</summary>
  <ol><li><b>Set up custom firewall</b>: find your current IP.  You can use this web page: https://whatismyipaddress.com/ to find your current IP.
Login to https://biohpc.cornell.edu/Default.aspx, under your user ID, click "Custom firewall". Use pulldown to locate your server, and enter your IP.</li>
 <li><b>Set up tunnel</b>:
Run this command on Mac Terminal or Windows Command Prompt. Replace USERNAME with your BioHPC user ID, cbsuXXXX with your assigned server, and 8009 with your mapped port (in both places):
<code>ssh -N -L 8009:localhost:8009 USERNAME@cbsuXXXXX.biohpc.cornell.edu</code>
 </li>
    <li>
After you run this command, it may look like the terminal hangs there. That is ok. Keep this terminal Wwndow open.
</li>
<li>
  Visit the jupyter notebook URL above, but replace the servername with the word 'localhost'. The URL should look something like: http://localhost:8009/?token=this1is2a3very4long5token6string.
</li>
</ol>
 </details>

### 6. Now you should be in jupyter!
You should see your files in the left-side pane (assuming folder icon on upper-right is chosen). You can navigate to the scRNA-WS24/Lessons directory, and double-click on the Week4-scanpy.ipynb file to start working through the jupyter notebook. 

### 7. Claim files when finished
Because it is running in docker, the files created by the notebook will not belong to you on the server. From the server, you can execute the command `docker1 claim` to 'claim' the files and fix the permissions.

