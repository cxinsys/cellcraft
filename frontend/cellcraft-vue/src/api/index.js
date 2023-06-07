import axios from "axios";
import { setInterceptors } from "@/api/common/interceptors";

// axios.defaults.baseURL = 'http://127.0.0.1:8000' // 서버주소
// axios.defaults.headers.post['Content-Type'] = 'application/json;charset=utf-8'
// axios.defaults.headers.post['Access-Control-Allow-Origin'] = 'http://127.0.0.1:8000'
// axios.defaults.withCredentials = true

function createInstance() {
  const instance = axios.create({
    baseURL: "http://127.0.0.1:8000",
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

function exportData(data) {
  return instance.post("/routes/workflow/compile", data);
}

function getResult(filename) {
  return instance.post("/routes/workflow/result", filename);
}

function uploadForm(formData) {
  // FormData의 value 확인
  // for (let value of formData.values()) {
  //   console.log(value)
  // }
  return instance.post("/routes/files/upload", formData, {
    headers: {
      "Content-Type": "multipart/form-data",
    },
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

export {
  registerUser,
  loginUser,
  getUser,
  exportData,
  getResult,
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
};
