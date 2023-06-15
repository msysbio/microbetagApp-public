
## Ports 

I changed the mysql port of the vm to 4545 so I can have the 3306 for the container


## My database

We need to build the whole database in the vm or wherever (even on a laptop) and then export it to a `dbtest_data.sql` file. 
Here is how: 
```bash
 mysqldump -u [user name] –p [password] [database_name] > [dumpfilename.sql]

```

## To build the SSL keys 

Follow the instructions [here](https://www.baeldung.com/openssl-self-signed-cert)
and remove the pass phrase as shown [here](https://help.cloud66.com/docs/security/remove-passphrase).

> `openssl rsa -in [original.key] -out [new.key]`


## Build a Docker app with HTTPS 

### Option A: The `my_image` Docker image

You may have a look at the `Dockerfile.ssl_example` to see how one can come with a Docker image including the SSL files. 

```
docker run -p 443:443 my_image
```


### Option B: The `docker-compose.yml` alternative


To run this case:
```
docker-compose up
```

This option allows us to fire up several apps at the same time and this our database too. 





## Links

[Setting up and Running a MySQL Container](https://www.baeldung.com/ops/docker-mysql-container)



