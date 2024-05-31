# Step 1: Use nginx to serve the built files
FROM nginx:alpine
COPY dist /usr/share/nginx/html
COPY nginx.conf /etc/nginx/nginx.conf
EXPOSE 8002
CMD ["nginx", "-g", "daemon off;"]
