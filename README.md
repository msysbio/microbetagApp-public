
## Env

I made a new file under  `/etc/docker/` called `daemon.json`:

```bash
u0156635@gbw-s-msysbio01:/etc/docker$ more daemon.json 
{
  "data-root": "/data/docker"

}
```


## Ports 

I changed the mysql port of the vm to 4545 so I can have the 3306 for the container


## My database

We need to build the whole database in the vm or wherever (even on a laptop) and then export it to a `dbtest_data.sql` file. 
Here is how: 
```bash
 mysqldump -u [user name] –p [password] [database_name] > [dumpfilename.sql]

```

### To build the SSL keys 

Follow the instructions [here](https://www.baeldung.com/openssl-self-signed-cert)
and remove the pass phrase as shown [here](https://help.cloud66.com/docs/security/remove-passphrase).

> `openssl rsa -in [original.key] -out [new.key]`


## Initiate the Docker app with HTTPS 

Init a tmux session: 

```bash
tmux new -s dockerApp
```

From within the tmux session, fire the `docker-compose` command:

```bash
docker-compose up
```

This initiates 3 Docker containers at once: 
- a database 
- an ngninx 
- the actual app 

You should have now something like this:

![tmux docker](figs/init-app.png)


**Attention!**

You need to **detach** from the tmux session; **not** exit! 
If you exit, the `docer-compose` command will be shut and your app will go down. 
To detach, you need to press:
`ctrl+b` and then `d`.

You can attach to your running session any time by runnint 
```bash
tmux attach -t dockerApp
```






## Links

- [Setting up and Running a MySQL Container](https://www.baeldung.com/ops/docker-mysql-container)
- [Creating a Self-Signed Certificate With OpenSS](https://www.baeldung.com/openssl-self-signed-cert)

Challenges:
[`argparse` in nginx](https://github.com/benoitc/gunicorn/issues/1867)

