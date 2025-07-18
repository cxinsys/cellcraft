import axios from "axios";
import { setInterceptors } from "@/api/common/interceptors";

// axios.defaults.baseURL = 'http://127.0.0.1:8000' // 서버주소
// axios.defaults.headers.post['Content-Type'] = 'application/json;charset=utf-8'
// axios.defaults.headers.post['Access-Control-Allow-Origin'] = 'http://127.0.0.1:8000'
// axios.defaults.withCredentials = true

function createInstance() {
  const instance = axios.create({
    // baseURL: "http://0.0.0.0:8000",
    baseURL: process.env.VUE_APP_BASE_URL,
  });
  return setInterceptors(instance);
}

const instance = createInstance();

function registerUser(userData) {
  return instance.post("/routes/auth/register", userData);
}

function loginUser(userData) {
  const formData = new URLSearchParams();
  formData.append('username', userData.username);
  formData.append('password', userData.password);
  
  return instance.post("/routes/auth/login/access-token", formData.toString(), {
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
    },
  });
}

function getUser() {
  return instance.get("/routes/auth/me");
}

function getFilteredUsers(conditions) {
  return instance.get("/routes/admin/users", { params: conditions });
}

function getFilteredFiles(conditions) {
  return instance.get("/routes/admin/files", { params: conditions });
}

function getFilteredWorkflows(conditions) {
  return instance.get("/routes/admin/workflows", { params: conditions });
}

function getFilteredTasks(conditions) {
  return instance.get("/routes/admin/tasks", { params: conditions });
}

function getFilteredPlugins(conditions) {
  return instance.get("/routes/admin/plugins", { params: conditions });
}

function getUsersCount() {
  return instance.get("/routes/admin/users_count");
}

function getFilesCount() {
  return instance.get("/routes/admin/files_count");
}

function getWorkflowsCount() {
  return instance.get("/routes/admin/workflows_count");
}

function getTasksCount() {
  return instance.get("/routes/admin/tasks_count");
}

function getPluginsCount() {
  return instance.get("/routes/admin/plugins_count");
}

function getSystemResources() {
  return instance.get("/routes/admin/system/stats");
}

function exportData(data) {
  return instance.post("/routes/workflow/compile", data);
}

function getResult(WorkflowResult) {
  return instance.post("/routes/workflow/result", WorkflowResult, {
    responseType: "blob", // 서버로부터 받은 데이터를 blob 형태로 처리
  });
}

function getResults(WorkflowResult) {
  return instance.post("/routes/workflow/results", WorkflowResult);
}

function runVisualization(WorkflowUpdate) {
  return instance.post("/routes/workflow/visualization", WorkflowUpdate);
}

function getVisualizationResult(WorkflowResult) {
  return instance.post("/routes/workflow/visualization/result", WorkflowResult);
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
  return instance.post("/routes/files/find", fileInfo);
}

function findFolder(folder) {
  return instance.post("/routes/files/folder", folder);
}

function deleteFile(file) {
  return instance.post("/routes/files/delete", file);
}

function convertFile(file) {
  return instance.post("/routes/files/convert", file);
}

function checkOptions() {
  return instance.get("/routes/files/setup/check");
}

function getOptions(file) {
  return instance.get(`/routes/files/setup/${file}`);
}

function checkConvert(file) {
  return instance.get(`/routes/files/check/${file}`);
}

function getColumns(fileInfo) {
  return instance.post("/routes/files/columns", fileInfo);
}

function getClusters(fileInfo) {
  return instance.post("/routes/files/clusters", fileInfo);
}

function setupAlgorithm(options) {
  return instance.post("/routes/files/setup", options);
}

function getFileData(filename) {
  return instance.get(`/routes/files/data/${filename}`);
}

function getWorkflows() {
  return instance.get("/routes/workflow/me");
}

function findWorkflow(workflowInfo) {
  return instance.post("/routes/workflow/find", workflowInfo);
}

function saveWorkflow(workflow) {
  return instance.post("/routes/workflow/save", workflow);
}

function deleteWorkflow(workflow) {
  return instance.post("/routes/workflow/delete", workflow);
}

function saveWorkflowNodeFile(workflowNodeFileInfo) {
  return instance.post("/routes/workflow/node/save", workflowNodeFileInfo);
}

function deleteWorkflowNodeFile(workflowNodeFileInfo) {
  return instance.post("/routes/workflow/node/delete", workflowNodeFileInfo);
}

function readWorkflowNodeFile(workflowNodeFileInfo) {
  return instance.post("/routes/workflow/node/read", workflowNodeFileInfo);
}

function taskMonitoring(taskId) {
  return instance.get(`/routes/task/info/${taskId}`);
}

function userTaskMonitoring() {
  return instance.get("/routes/task/monitoring");
}

function revokeTask(taskId) {
  return instance.delete(`/routes/task/revoke/${taskId}`);
}

function deleteTask(taskId) {
  return instance.delete(`/routes/task/delete/${taskId}`);
}

function getTaskLogs(taskId) {
  return instance.get(`/routes/task/logs/${taskId}`);
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

function getTutorialFileDownload(filename) {
  return instance.get(`/routes/files/tutorials/${filename}`, {
    responseType: "blob", // 서버로부터 받은 데이터를 blob 형태로 처리
  });
}

function getResultFileOne(filename) {
  return instance.get(`/routes/files/result/${filename}`);
}

function validationPlugin(plugin, rules, drawflow) {
  const data = {
    plugin,
    rules,
    drawflow,
  };

  return instance.post("/routes/plugin/validation", data);
}

function uploadPluginMetadata(pluginCreate) {
  return instance.post("/routes/plugin/upload", pluginCreate);
}

function syncPluginData(pluginName) {
  return instance.post(`/routes/plugin/update/${pluginName}`);
}

function uploadPluginScripts(formData) {
  return instance.post("/routes/plugin/upload_scripts", formData, {
    headers: {
      "Content-Type": "multipart/form-data",
    },
  });
}

function uploadPluginPackage(formData) {
  return instance.post("/routes/plugin/upload_package", formData, {
    headers: {
      "Content-Type": "multipart/form-data",
    },
  });
}

function uploadTextDependencies(formData) {
  return instance.post("/routes/plugin/upload_text_dependencies", formData);
}

function getPlugins() {
  return instance.get("/routes/plugin/list");
}

function associatePlugin(plugin_id) {
  return instance.post("/routes/plugin/associate", { plugin_id : plugin_id });
}

function dissociatePlugin(plugin_id) {
  return instance.post("/routes/plugin/dissociate", { plugin_id : plugin_id });
}

function getPluginTemplate(plugin_id) {
  return instance.get(`/routes/plugin/template/${plugin_id}`);
}

function getPluginFile(file_info) {
  return instance.get(`/routes/plugin/file/${file_info.plugin_name}/${file_info.file_name}`, {
    responseType: "blob", // 서버로부터 받은 데이터를 blob 형태로 처리
  });
}

function getPluginPackageList(plugin_name) {
  return instance.get(`/routes/plugin/package/${plugin_name}`);
}

function getPluginReferenceFolders(plugin_name) {
  return instance.get(`/routes/plugin/reference_folders/${plugin_name}`);
}

function getDataTableFile(vgt_info) {
  return instance.post("/routes/datatable/load_data", vgt_info);
}

function getPluginInfo(name) {
  return instance.get(`/routes/plugin/info/${name}`);
}

function buildPlugin(plugin_name) {
  return instance.post(`/routes/plugin/build/${plugin_name}`);
}

function buildPluginDocker(plugin_name, use_gpu = false) {
  return instance.post(`/routes/plugin/build_docker/${plugin_name}`, { use_gpu });
}

function checkPluginImage(plugin_name) {
  return instance.get(`/routes/plugin/check_image/${plugin_name}`);
}

// 새로운 플러그인 빌드 관련 API 함수들
function getBuildStatus(task_id) {
  return instance.get(`/routes/plugin/build/status/${task_id}`);
}

function getBuildTasks(plugin_name = null) {
  const params = plugin_name ? { plugin_name } : {};
  return instance.get('/routes/plugin/build/tasks', { params });
}

function cancelBuildTask(task_id) {
  return instance.post(`/routes/plugin/build/cancel/${task_id}`);
}

function getBuildLogs(plugin_name) {
  return instance.get(`/routes/plugin/build/logs/${plugin_name}`);
}

function createTaskEventSource(taskId, callbacks = {}) {
  const eventSource = new EventSource(`${process.env.VUE_APP_BASE_URL}/routes/task/info/${taskId}`);
  
  eventSource.onmessage = (event) => {
    // 기본 메시지 핸들러
    if (callbacks.onMessage) {
      callbacks.onMessage(event);
    }
    
    // 작업 완료시 이벤트소스 자동 종료
    if (event.data === "SUCCESS" || event.data === "FAILURE" || event.data === "REVOKED") {
      if (callbacks.onComplete) {
        callbacks.onComplete(event.data);
      }
      eventSource.close();
    }
  };

  // 에러 핸들링
  eventSource.onerror = (error) => {
    if (callbacks.onError) {
      callbacks.onError(error);
    }
    eventSource.close();
  };

  return eventSource;
}

// Admin APIs
function updateUserAdmin(userId, userData) {
  return instance.put(`/routes/admin/users/${userId}`, userData);
}

function deleteUserAdmin(userId) {
  return instance.delete(`/routes/admin/users/${userId}`);
}

function deleteFileAdmin(fileId) {
  return instance.delete(`/routes/admin/files/${fileId}`);
}

function deleteWorkflowAdmin(workflowId) {
  return instance.delete(`/routes/admin/workflows/${workflowId}`);
}

export {
  registerUser,
  loginUser,
  getUser,
  getUsersCount,
  getFilesCount,
  getWorkflowsCount,
  getTasksCount,
  getPluginsCount,
  getFilteredUsers,
  getFilteredFiles,
  getFilteredWorkflows,
  getFilteredTasks,
  getFilteredPlugins,
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
  deleteTask,
  getTaskLogs,
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
  validationPlugin,
  uploadPluginMetadata,
  syncPluginData,
  uploadPluginScripts,
  uploadPluginPackage,
  uploadTextDependencies,
  getPlugins,
  associatePlugin,
  dissociatePlugin,
  getPluginTemplate,
  getDataTableFile,
  saveWorkflowNodeFile,
  deleteWorkflowNodeFile,
  readWorkflowNodeFile,
  getFileData,
  getSystemResources,
  getPluginFile,
  getPluginPackageList,
  getPluginInfo,
  getVisualizationResult,
  runVisualization,
  createTaskEventSource,
  getPluginReferenceFolders,
  getTutorialFileDownload,
  updateUserAdmin,
  deleteUserAdmin,
  deleteFileAdmin,
  deleteWorkflowAdmin,
  buildPlugin,
  buildPluginDocker,
  checkPluginImage,
  getBuildStatus,
  getBuildTasks,
  cancelBuildTask,
  getBuildLogs,
};
