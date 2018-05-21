#! /bin/bash
sudo apt-get update
sudo apt-get install git libuv1-dev libssl-dev gcc g++ cmake make -y
git clone https://github.com/uWebSockets/uWebSockets 
cd uWebSockets
git checkout e94b6e1
mkdir build
cd build
cmake ..
make 
make install
cd ../..
ln -s /usr/lib64/libuWS.so /usr/lib/libuWS.so
rm -r uWebSockets
