#Installation

mkdir ~/Desktop/docker
cd ~/Desktop/docker
sudo apt-get install \
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
apt-cache madison docker-ce
sudo apt-get install docker-ce=5:20.10.7~3-0~ubuntu-xenial docker-ce-cli=5:20.10.7~3-0~ubuntu-xenial containerd.io
sudo docker run hello-world

#Check status
sudo service docker status

docker.service - Docker Application Container Engine
   Loaded: loaded (/lib/systemd/system/docker.service; enabled; vendor preset: enabled)
   Active: active (running) since Tue 2021-11-02 11:32:10 EDT; 1min 56s ago
     Docs: https://docs.docker.com
 Main PID: 29573 (dockerd)
    Tasks: 72
   Memory: 93.6M
      CPU: 607ms
   CGroup: /system.slice/docker.service
           └─29573 /usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock
           
#Install shinyproxy         
git clone https://github.com/openanalytics/shinyproxy.git
sudo apt install maven
cd shinyproxy
mvn -U clean install

#Create directory and place app files in it
mkdir seurat_docker
mkdir seurat_docker/app

#within app dir, place server.R, ui.R, functions.R and data folder
#withing seurat_docker dir, place Dockerfile

#create image
cd seurat_docker
sudo docker build -t seurat_docker .

cd ngs_docker
sudo docker build -t ngs_docker .

#view docker imager
sudo docker images

docker version

#Run shiny-proxy
java -jar shinyproxy/shinyproxy-2.5.0.jar

#Upload to dockerhub
#Create repository ngs_viewer

docker tag ngs_docker:latest apoorvababu/ngs_viewer:latest
docker push apoorvababu/ngs_viewer:latest

#Create repository seurat_viewer

docker tag seurat_docker:latest apoorvababu/seurat_viewer:latest
docker push apoorvababu/seurat_viewer:latest
