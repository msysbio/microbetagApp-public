FROM nginx:1.25

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d
RUN sed -E -i '17i \ \ \ \ client_max_body_size 20M;' /etc/nginx/nginx.conf  
    # && \
    # sed -E -i '' /etc/nginx/nginx.conf

RUN apt-get update && apt-get install -y apt-transport-https
RUN apt-get install -y vim