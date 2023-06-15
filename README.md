

# To build the SSL keys 

Follow the instructions [here](https://www.baeldung.com/openssl-self-signed-cert)
and remove the pass phrase as shown [here](https://help.cloud66.com/docs/security/remove-passphrase).

> `openssl rsa -in [original.key] -out [new.key]`


# Build a Docker app with HTTPS 

## Option A: The `my_image` Docker image

You may have a look at the `Dockerfile.ssl_example` to see how one can come with a Docker image including the SSL files. 

```
docker run -p 443:443 my_image
```


## Option B: The `docker-compose.yml` alternative


To run this case:
```
docker-compose up
```

This option allows us to fire up several apps at the same time and this our database too. 



