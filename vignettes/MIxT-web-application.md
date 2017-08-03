# MIxT web application
The MIxT web application is designed for exploring the results from the MIxT
analysis comparing transcriptional profiles from two matched tissues across
individuals. 

For a more detailed description of the ideas and design of the MIxT web
application please refer to
["Building Applications For Interactive Data Exploration In Systems Biology" by Fjukstad et al.](biorxiv.org/content/early/2017/05/24/141630) 
and the source code [here](https://github.com/fjukstad/mixt). 

In summary the application consists of three components or services: a
compute service that provides analyses of transcriptional profiles, a database
service for retrieving metadata for e.g. genes, and a web server hosting the
application itself. The web server interface with both the compute and database
services to retrieve data in the web application. 

If you  want to run the web application you must run all services. We have
bundled them all in [Docker](http://docker.com) containers that make deploying
the application a simple task: 

```
$ git clone github.com/fjukstad/mixt
$ cd mixt
$ docker-compose up
```

The application should now run on [localhost:8000](http://localhost:8000). The
above workflow requires that you have both `git` and [Docker](http://docker.com)
installed. 

# New Data 
The [mixtApp](https://github.com/vdumeaux/mixtApp) R package provides data and
analyses of transcriptional profiles for the web application. Instead of
building a web application where we use pre-computed results we run the analyses
on demand. 

If you want to modify the data being used in the web application you'll have to
rebuild the mixtApp package yourself, and start the different services
manually.

First clone down the  [mixtApp](https://github.com/vdumeaux/mixtApp)
repository: 

```
$ git clone https://github.com/vdumeaux/mixtApp.git
```

Then `cd` into the `mixtApp` directory and add your data to the `data/` folder
formatted as described earlier. Next up is building the compute service
container with your new mixtApp package.

Still in the `mixtApp/` directory, run: 

```
docker build -t compute-service .
```

which build the container for you. Now you can start it up to accept requests on
port `8787` by running 

```
docker run -p --name=compute-service -t compute-service
```

and it should appear with the `docker ps` command: 

```
 docker ps
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS               NAMES
08c9889b705e        compute-service     "/bin/sh -c 'go ru..."   4 seconds ago       Up 4 seconds        80/tcp              compute-service
```

The compute service is now running, so next up is starting the web application
container. This is lucily one liner: 

```
docker run -p 8000:80 --link compute-service -e COMPUTE_SERVICE=compute-service:80 --name=mixt -t fjukstad/mixt
```

That's it!  You can now visit the application running on
[localhost:8000](http://localhost:8000). 

If you need more details on the docker commands you can have a look
at[docs.docker.com](https://docs.docker.com).


