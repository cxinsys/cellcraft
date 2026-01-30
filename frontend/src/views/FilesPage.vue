<template>
  <div class="layout" @click="ClickOut">
    <section class="folder">
      <div class="folder__header">
        <p class="folder__title left">Navigator</p>
        <div class="folder__createBtn right">
          <!-- <img class="folder__createBtn--icon" src="@/assets/add-folder.png" />
          <img class="folder__createBtn--icon" src="@/assets/refresh.png" /> -->
        </div>
      </div>
      <ul class="folder__list">
        <li class="folder__item" v-for="(folder, idx) in folders_list" :key="idx"
          v-bind:class="{ toggleFolder: toggleFolder === idx }" @click="folderClick(idx, folder[0])">
          <div class="folder__item--col">
            <img class="folder__item--icon" src="@/assets/arrow-bottom.png" v-if="toggleFolder === idx" />
            <img class="folder__item--icon" src="@/assets/arrow-right.png" v-else />
            <img class="folder__item--icon large" src="@/assets/open-folder.png" v-if="toggleFolder === idx" />
            <img class="folder__item--icon large" src="@/assets/folder.png" v-else />
          </div>
          <p class="folder__name">{{ folder[0] }}</p>
        </li>
      </ul>
    </section>
    <main class="files">
      <div class="files__header">
        <div class="header__column left">
          <p class="files__folder">{{ currentFolder }}</p>
        </div>
        <!-- <input
          class="files__search"
          type="text"
          name="search"
          placeholder="search anything..."
        /> -->
        <div class="header__column right">
          <label class="files__button">
            <img class="files__button--icon" src="@/assets/upload-file-black.png" />
            <h1>Upload File</h1>
            <input class="files__input" type="file" name="file" ref="selectFile" @change.prevent="uploadFile" />
          </label>
          <div class="progress__box" v-if="uploadPercentage > 0">
            <progress :value="uploadPercentage" max="100"></progress>
            <span>{{ uploadPercentage }}%</span>
          </div>
        </div>
      </div>

      <table class="files__table">
        <thead>
          <tr>
            <th>Name</th>
            <th>Date</th>
            <th>Type</th>
            <th>Size</th>
          </tr>
        </thead>
        <tbody>
          <tr class="files__item" v-for="(file, idx) in files_list" :key="idx" @contextmenu.prevent
            @click.right="RMouseClick($event, file.file_name, idx)" v-bind:class="{ select: R_Mouse_isActive }">
            <td>{{ extractFileName(file.file_name) }}</td>
            <td>{{ cutDateFromISO(file.created_at) }}</td>
            <td>{{ extractExtension(file.file_name) }}</td>
            <td>{{ formatBytes(file.file_size) }}</td>
          </tr>
        </tbody>
      </table>
      <ul class="toggle__menu" v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }">
        <!-- <li>view</li>
        <li>plot</li>
        <li>rename</li> -->
        <li @click="removeFile">Delete</li>
      </ul>
    </main>
  </div>
</template>

<script>
import { uploadForm, getFiles, findFolder, deleteFile } from "@/api/index";
import { validateFileExtension } from "@/utils/validation";
import { generateUploadFileName } from "@/utils/filename";
import { calculateContextMenuPosition } from "@/utils/positionCalculator";
import { formatBytes, cutDateFromISO, extractFileName, extractExtension } from "@/utils/formatters";

export default {
  props: {
    doUploadFile: {
      type: String,
      default: "",
    },
  },

  data() {
    return {
      folders_list: [],
      files_list: [],
      R_Mouse_isActive: false,
      Clickout_isActive: false,
      toggleFolder: null,
      currentFolder: "data",
      xPosition: 0,
      yPosition: 0,
      selectFile: null,
      file_name: null,
      list_idx: null,
      uploadPercentage: 0,
      isDeletingFile: false, // Track file deletion state to prevent concurrent deletions
    };
  },

  methods: {
    // Formatter methods for template
    formatBytes,
    cutDateFromISO,
    extractFileName,
    extractExtension,

    RMouseClick(event, file_name, idx) {
      this.R_Mouse_isActive = false;
      const { x, y } = calculateContextMenuPosition(event.clientX, event.clientY);
      this.xPosition = x;
      this.yPosition = y;
      this.R_Mouse_isActive = true;
      this.file_name = file_name;
      this.list_idx = idx;
    },
    ClickOut() {
      this.R_Mouse_isActive = false;
    },
    folderClick(idx, folderName) {
      if (idx === this.toggleFolder) {
        this.toggleFolder = null;
      } else {
        this.toggleFolder = idx;
      }
      this.currentFolder = folderName;
    },
    async uploadFile() {
      if (this.$refs.selectFile.files.length > 0) {
        const file = this.$refs.selectFile.files[0];

        // Validate file extension using utility function
        const validation = validateFileExtension(file.name);
        if (!validation.isValid) {
          alert(validation.message);
          return;
        }

        // Generate upload filename using utility function
        const newFileName = generateUploadFileName(this.currentFolder, file.name);
        this.selectFile = new File([file], newFileName);

        const form = new FormData();
        form.append("files", this.selectFile);

        // 파일 업로드 진행률을 추적하기 위한 콜백
        const onUploadProgress = (progressEvent) => {
          this.uploadPercentage = parseInt(
            Math.round((progressEvent.loaded * 100) / progressEvent.total)
          );
        };

        try {
          await uploadForm(form, onUploadProgress);
          this.uploadPercentage = 0; // 업로드 완료 후 초기화

          // 업로드 성공 후 파일 목록 새로고침
          const folderList = await findFolder({
            folder_name: this.currentFolder,
          });
          this.files_list = folderList.data;
        } catch (error) {
          console.error('File upload failed:', error);

          // 업로드 실패 시 진행률 리셋
          this.uploadPercentage = 0;

          // 사용자에게 에러 메시지 표시
          let errorMessage = 'File upload failed. Please try again.';

          if (error.response) {
            // 백엔드 에러 응답이 있는 경우
            if (error.response.data && error.response.data.detail) {
              errorMessage = error.response.data.detail;
            } else if (error.response.status === 413) {
              errorMessage = 'File is too large. Please choose a smaller file.';
            } else if (error.response.status === 400) {
              errorMessage = 'Invalid file. Please check the file format and try again.';
            } else if (error.response.status === 500) {
              errorMessage = 'Server error occurred. Please try again later.';
            }
          } else if (error.request) {
            errorMessage = 'Network error. Please check your connection.';
          }

          alert(errorMessage);
        }
      }
    },
    async removeFile() {
      // Prevent concurrent deletions
      if (this.isDeletingFile) {
        return;
      }

      // Show confirmation dialog BEFORE any state changes
      if (
        !confirm(
          "Are you sure you want to delete this file? This action cannot be undone."
        )
      ) {
        return; // User cancelled
      }

      // Set loading state
      this.isDeletingFile = true;

      try {
        const file = {
          file_name: this.file_name,
        };

        // Call API and wait for response
        await deleteFile(file);

        // Only remove from UI after successful deletion
        this.files_list.splice(this.list_idx, 1);

      } catch (error) {
        console.error('File deletion failed:', error);

        // User-friendly error messages
        let errorMessage = 'Failed to delete file. Please try again.';
        if (error.response?.data?.detail) {
          errorMessage = error.response.data.detail;
        } else if (error.response?.status === 404) {
          errorMessage = 'File not found. It may have been already deleted.';
        } else if (error.request) {
          errorMessage = 'Network error. Please check your connection.';
        }

        alert(errorMessage);
      } finally {
        // Reset loading state
        this.isDeletingFile = false;
      }
    },
  },

  async mounted() {
    if (this.$route.query.doUploadFile) {
      this.$refs.selectFile.click();
    }
    try {
      const fileList = await getFiles();
      this.folders_list = fileList.data;
      const folderList = await findFolder({
        folder_name: this.currentFolder,
      });
      this.files_list = folderList.data;
    } catch (error) {
      console.error(error);
    }
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  position: relative;
  overflow: auto;
}

.folder {
  width: 18rem;
  height: 100%;
  border-right: 1px solid #aeaeae;
  background-color: rgb(201, 202, 203);
  /* width: 22rem;
  height: 100%;
  background: #cfcfcf;
  border-right: 1px solid #afafaf; */
}

.folder__header {
  width: 90%;
  height: 10%;
  margin: auto;
  display: flex;
  align-items: center;
  position: relative;
  color: rgba(0, 0, 0, 0.8);
}

.folder__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.7rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
}

.folder__createBtn {
  width: 6rem;
  height: 2rem;
  vertical-align: middle;
}

.folder__createBtn--icon {
  width: 2rem;
  height: 2rem;
  object-fit: cover;
  margin: 0 0.5rem;
  opacity: 0.8;
}

.folder__list {
  width: 85%;
  height: 90%;
  margin: auto;
}

.folder__item {
  width: 100%;
  height: 5%;
  /* margin: 0 0rem; */
  border-radius: 0.5rem;
  display: flex;
  cursor: pointer;
}

.folder__item:hover {
  background: rgb(176, 177, 178);
}

.toggleFolder {
  background: rgb(176, 177, 178);
}

.folder__item--col {
  width: 5rem;
  height: 100%;
  display: flex;
  align-items: center;
}

.folder__item--icon {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
}

.large {
  width: 1.5rem;
  height: 1.5rem;
}

.folder__name {
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
}

.files {
  width: 75%;
  height: 100%;
  background: #ffffff;
}

.files__header {
  width: 100%;
  height: 10%;
  display: flex;
  align-items: center;
  justify-content: center;
  position: relative;
}

.header__column {
  width: 50%;
  height: 100%;
  display: flex;
  align-items: center;
}

.right {
  position: absolute;
  right: 0;
}

.left {
  position: absolute;
  left: 0;
}

.files__folder {
  margin-left: 2rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.8rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
  /* text-transform: capitalize; */
}

.files__button {
  width: 9rem;
  height: 2rem;
  padding: 0.2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  /* background: #ffffff; */
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.files__button:hover {
  cursor: pointer;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.files__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
}

.files__input {
  display: none;
}

.files__search {
  width: 300px;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 2rem;
  outline-style: none;
  background: #f7f7f7;
}

.files__search:focus {
  border: 1px solid #bcbcbc;
}

.files__table {
  width: 95%;
  height: auto;
  margin: auto;
  border-collapse: collapse;
}

.files__table thead {
  height: 3rem;
  border-bottom: 1px solid #e1e1e1;
}

.files__table th {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.1rem;
  line-height: 1.1rem;

  text-align: left;
  vertical-align: middle;
  padding: 0 1rem;

  color: rgba(0, 0, 0, 0.5);
}

.files__table td {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  vertical-align: middle;
  padding: 1rem;
}

.files__item:hover {
  cursor: pointer;
  background: rgb(204, 218, 245);
}

.files__item td {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.files__header__img {
  margin-left: 5px;
  width: 15px;
}

.font-setting {
  font-size: 1vw;
  font-weight: bold;
  color: rgb(70, 70, 70);
}

.files__contents {
  width: 100vw;
  height: 75vh;

  margin-top: 20px;
}

.files__contents__item {
  display: flex;
  align-items: center;

  width: 100%;
  height: 40px;
}

/* .files__contents__item.select{
    background-color: #242F9B;
    color: white;
  } */

.files__contents__item:hover {
  background-color: #242f9b;
  color: white;
}

.toggle__menu {
  display: none;
  position: absolute;
  width: 200px;
  margin: 0;
  padding: 0;
  background: #ffffff;
  border-radius: 5px;
  list-style: none;
  box-shadow: 0 15px 35px rgba(50, 50, 90, 0.1), 0 5px 15px rgba(0, 0, 0, 0.07);
  overflow: hidden;
  z-index: 999999;
  text-transform: capitalize;
}

.toggle__menu.open {
  display: block;
  opacity: 1;
  position: absolute;
}

.toggle__menu>li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}

.toggle__menu>li:hover {
  background: #e5e5e5;
}

.progress__box {
  display: flex;
  align-items: center;
  justify-content: center;
  width: 15rem;
  height: 1rem;
}

/* progress bar */
progress {
  width: 100%;
  /* 전체 너비를 차지하도록 설정 */
  height: 100%;
  /* 높이 설정 */
  background-color: #eee;
  /* 배경색 설정 */
  border-radius: 10px;
  /* 모서리 둥글게 처리 */
  margin-right: 0.5rem;
}

progress::-webkit-progress-bar {
  background-color: #eee;
  /* 크롬, 사파리 등 WebKit 기반 브라우저에서의 배경색 */
}

progress::-webkit-progress-value {
  background-color: #4caf50;
  /* 크롬, 사파리 등 WebKit 기반 브라우저에서의 진행률 색상 */
}

progress::-moz-progress-bar {
  background-color: #4caf50;
  /* 파이어폭스에서의 진행률 색상 */
}
</style>
