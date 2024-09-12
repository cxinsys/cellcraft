# Step 1: Build the Vue.js application
FROM node:14 AS build-stage

WORKDIR /app

# Install dependencies
COPY package*.json ./
RUN npm install

# Copy source files and build
COPY . .
RUN npm run build

# Step 2: Use Nginx to serve the built files
FROM nginx:alpine

# Copy built files from build-stage
COPY --from=build-stage /app/dist /usr/share/nginx/html

# Copy custom nginx configuration
COPY nginx.conf /etc/nginx/nginx.conf

EXPOSE 8080

CMD ["nginx", "-g", "daemon off;"]
