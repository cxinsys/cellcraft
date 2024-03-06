FROM node:alpine as build-stage

WORKDIR /app

COPY package*.json ./

RUN npm install

COPY . .

RUN npm run build

FROM nginx:alpine

COPY --from=build-stage /app/dist /usr/share/nginx/html

# Path: frontend/nginx/default.conf
COPY nginx/default.conf /etc/nginx/conf.d/default.conf