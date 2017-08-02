# MIxT web application
The MIxT web application is designed for exploring the results from the MIxT analysis comparing transcriptional profiles from two matched tissues across individuals. 


##  mixtApp R Package (mixtApp)
R package that provides the extra functions used in the compute service for
the MIxT web application. 

### Install 
Using R 
```
> devtoools::install_github("vdumeaux/mixtApp")
```
or with the shell 
```
$ git clone https://github.com/vdumeaux/mixtApp.git
$ R CMD INSTALL mixtApp
```
if you plan on using it with your data later. 

### Data ?
If you want to build the mixt web application all data should be formatted as 
decribed earlier and saved in `data/` folder of the mixtApp package. 

You can modify `data-raw/datasets.R` to retrieve data and place it in the
`data/` folder in the R package. Then reinstall the package.

```
$ git clone git@github.com:vdumeaux/mixtApp.git 
$ cd mixtApp
# modify the data-raw/datasets.R file to load your data. 
$ R -f data-raw/datasets.R
$ R CMD INSTALL .
```

## Backend data and analysis server for MIxT 
The presentation (web app) and data analysis (compute-backend) are separated into two services/processes/docker containers


First, install [Docker](http://docker.com) 


### Build and run the compute service

This the is the compute service in the MIxT web application. It is simply the
compute service in [Kvik](https://github.com/fjukstad/kvik) with the
[mixtApp](http://github.com/vdumeaux/mixtApp) package installed. 

You can clone the repository, 
```
git clone git@github.com:fjukstad/mixt-compute-service.git
```
To use mixtApp R package witht he new data, you need to modify the [Dockerfile]
(https://github.com/fjukstad/mixt-compute-service/blob/master/Dockerfile)
so it points to your new mixtApp Rpackage then build the Docker image, 
and run the container:

```
docker build -t mixt-compute-service .
docker run -p 8787:80 -t mixt-compute-service
```
Note: If there's already a server running you'll have to use `docker ps` to get
its name and `docker stop CONTAINTERNAME` to stop it.

### Build and run the web application

The compute service runs on port `:8787` and can be used by the
[MIxT web application](http://github.com/fjukstad/mixt). 

First install [go](https://golang.org/)
and create a new branch in the mixt directory

```
go get github.com/fjukstad/mixt
$ cd $GOPATH/src/github.com/fjukstad/mixt
git checkout -b myApp
```

in this branch modify the following files

then build the docker image and run the container
```

```
