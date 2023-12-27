import axios from "axios";
import { setInterceptors } from "@/api/common/interceptors";

// axios.defaults.baseURL = 'http://127.0.0.1:8000' // 서버주소
// axios.defaults.headers.post['Content-Type'] = 'application/json;charset=utf-8'
// axios.defaults.headers.post['Access-Control-Allow-Origin'] = 'http://127.0.0.1:8000'
// axios.defaults.withCredentials = true

function createInstance() {
  const instance = axios.create({
    // baseURL: "http://localhost:8082",
    baseURL: "http://cellcraft.ssu.ac.kr:8000",
    // baseURL: "http://backend:8000",
    // baseURL: "http://127.0.0.1:8000",
  });
  return setInterceptors(instance);
}

const instance = createInstance();

function registerUser(userData) {
  return instance.post("/routes/auth/register", userData);
}

function loginUser(userData) {
  return instance.post("/routes/auth/login/access-token", userData, {
    headers: {
      "Content-Type": "multipart/form-data",
    },
  });
}

function getUser() {
  return instance.get("/routes/auth/me");
}

function getFilteredUsers(conditions) {
  return instance.get("/routes/admin/users", { params: conditions });
}

function getUsersCount() {
  return instance.get("/routes/admin/users_count");
}

function exportData(data) {
  return instance.post("/routes/workflow/compile", data);
}

function getResult(filename) {
  return instance.post("/routes/workflow/result", filename);
}

function getResults() {
  return instance.post("/routes/workflow/results");
}

function uploadForm(formData, onUploadProgress) {
  // FormData의 value 확인
  // for (let value of formData.values()) {
  //   console.log(value)
  // }
  return instance.post("/routes/files/upload", formData, {
    headers: {
      "Content-Type": "multipart/form-data",
    },
    onUploadProgress,
  });
}

function getFiles() {
  return instance.get("/routes/files/me");
}

function findFile(fileInfo) {
  return instance.post("routes/files/find", fileInfo);
}

function findFolder(folder) {
  return instance.post("routes/files/folder", folder);
}

function deleteFile(file) {
  return instance.post("routes/files/delete", file);
}

function convertFile(file) {
  return instance.post("routes/files/convert", file);
}

function checkOptions() {
  return instance.get("routes/files/setup/check");
}

function getOptions(file) {
  return instance.get(`routes/files/setup/${file}`);
}

function checkConvert(file) {
  return instance.get(`routes/files/check/${file}`);
}

function getColumns(fileInfo) {
  return instance.post("/routes/files/columns", fileInfo);
}

function getClusters(fileInfo, annotation) {
  return instance.post("/routes/files/clusters", fileInfo, annotation);
}

function setupAlgorithm(options) {
  return instance.post("/routes/files/setup", options);
}

function getWorkflows() {
  return instance.get("routes/workflow/me");
}

function findWorkflow(workflowInfo) {
  return instance.post("routes/workflow/find", workflowInfo);
}

function saveWorkflow(workflow) {
  return instance.post("/routes/workflow/save", workflow);
}

function deleteWorkflow(workflow) {
  return instance.post("/routes/workflow/delete", workflow);
}

function taskMonitoring(taskId) {
  return instance.get(`/routes/workflow/task/${taskId}`);
}

function userTaskMonitoring() {
  return instance.get("/routes/workflow/monitoring");
}

function revokeTask(taskId) {
  return instance.delete(`/routes/workflow/revoke/${taskId}`);
}

function getHtml(filename) {
  return instance.get(`/routes/files/html/${filename}`);
}

function getResultFile(fileInfo) {
  return instance.post("/routes/files/result", fileInfo);
}

function getDownloadResult(filename) {
  return instance.get(`/routes/files/result/${filename}`, {
    responseType: "blob", // 서버로부터 받은 데이터를 blob 형태로 처리
  });
}

function getResultFileOne(filename) {
  return instance.get(`/routes/files/result/${filename}`);
}

export {
  registerUser,
  loginUser,
  getUser,
  getUsersCount,
  getFilteredUsers,
  exportData,
  getResult,
  getResults,
  uploadForm,
  getFiles,
  findFile,
  findFolder,
  getWorkflows,
  findWorkflow,
  saveWorkflow,
  deleteWorkflow,
  taskMonitoring,
  userTaskMonitoring,
  deleteFile,
  revokeTask,
  convertFile,
  getColumns,
  getClusters,
  setupAlgorithm,
  checkConvert,
  getOptions,
  checkOptions,
  getHtml,
  getResultFile,
  getDownloadResult,
  getResultFileOne,
};
