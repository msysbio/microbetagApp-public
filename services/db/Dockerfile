
FROM mysql:latest

RUN microdnf install -y vim
RUN microdnf install --assumeyes python3-pip
RUN pip3.9 install mysql-connector-python

RUN chown -R mysql:root /var/lib/mysql/ 

ARG MYSQL_DATABASE
ARG MYSQL_USER
ARG MYSQL_PASSWORD
ARG MYSQL_ROOT_PASSWORD

ENV MYSQL_DATABASE=$MYSQL_DATABASE
ENV MYSQL_USER=$MYSQL_USER
ENV MYSQL_PASSWORD=$MYSQL_PASSWORD
ENV MYSQL_ROOT_PASSWORD=$MYSQL_ROOT_PASSWORD


COPY init.sql /docker-entrypoint-initdb.d/

EXPOSE 3306
