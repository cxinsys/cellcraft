<template>
  <div class="layout">
    <form class="fileUpload-form" @submit.prevent="applyFile">
      <ul class="folder__list">
        <router-link class="cloud-row" target="_blank" to="/files">
          <div class="cloud-textbox">
            <p class="cloud-textbox__directory__desc">
              If you want to upload a new file, please click here.
            </p>
            <h1 class="cloud-textbox__directory">Files directory ></h1>
          </div>
        </router-link>
        <li class="folder__item" v-for="(folder, idx) in folders_list" :key="idx + 'O'"
          v-bind:class="{ toggleFolder: toggleFolder === idx }" @click="folderClick(idx, folder[0])">
          <div class="folder__item--col">
            <img class="folder__item--arrow" src="@/assets/arrow-right.png" v-if="toggleFolder === idx" />
            <img class="folder__item--arrow" src="@/assets/arrow-right.png" v-else />
            <img class="folder__item--icon" src="@/assets/open-folder.png" v-if="toggleFolder === idx" />
            <img class="folder__item--icon" src="@/assets/folder.png" v-else />
          </div>
          <p class="folder__name">{{ folder[0] }}</p>
        </li>
      </ul>
      <div class="fileUpload">
        <div class="form-row">
          <div class="form__selectFile">
            <ul class="form__fileList">
              <li class="file__item" v-for="(file, idx) in files_list" :key="idx + 'I'"
                v-bind:class="{ toggleFile: toggleFile === idx }" @click="fileClick(idx, file.file_name)">
                <div class="folder__item--col">
                  <img class="folder__item--arrow" src="@/assets/arrow-right.png" />
                  <img class="folder__item--icon" src="@/assets/file-icon.png" />
                </div>
                <p class="folder__name">{{ file.file_name }}</p>
              </li>
            </ul>
          </div>
        </div>
        <div class="form-row">
          <div class="form__name">Current File</div>
        </div>
        <div class="form-row">
          <ul class="form__info" v-if="selectFile">
            <li class="form__info--type">
              <img class="form__info--img" src="@/assets/file-icon.png" />
            </li>
            <li class="form__info--name">
              {{ selectFile.file_name }}&nbsp;&nbsp;&nbsp;{{
                selectFile.file_size | formatBytes
              }}
            </li>
          </ul>
          <ul class="form__info" v-else>
            <li class="form__info--blank">Please add data file</li>
          </ul>
          <label class="form__button--apply" v-bind:class="{ activate: apply }">
            {{ apply ? "Applied" : "Apply" }}
            <input class="form__input" type="submit" value="업로드" />
          </label>
          <div v-if="isLoading" class="loading-layout">
            <span> </span>
          </div>
        </div>
      </div>
    </form>
  </div>
</template>

<script>
import {
  getFiles,
  findFolder,
  findFile,
  convertFile,
  checkConvert,
} from "@/api/index";

export default {
  data() {
    return {
      node_name: "File",
      selectFile: null,
      is_upload: false,
      done_upload: false,
      getFile: false,
      toggleFolder: null,
      toggleFile: null,
      folders_list: [],
      files_list: [],
      recentFiles_list: [],
      apply: false,
      isLoading: false,
      nodeId: this.$route.query.node,
    };
  },
  methods: {
    previewFile() {
      if (this.$refs.selectFile.files.length > 0) {
        this.selectFile = this.$refs.selectFile.files[0];
      }
    },
    filterAndAddSuffix(inputString) {
      // Check if the inputString contains an underscore
      if (inputString.includes("_")) {
        // "_"로 구분된 문자열을 배열로 변환
        const segments = inputString.split("_");
        // 마지막 두 요소를 제외한 나머지를 합침
        const fileName = segments.slice(0, -2).join("_") + ".h5ad";
        return fileName;
      }
      // If no underscore found, return the original string
      return inputString;
    },
    async getFinder() {
      this.getFile = true;
      try {
        const fileList = await getFiles();
        this.folders_list = fileList.data;
      } catch (error) {
        console.error(error);
      }
    },
    async folderClick(idx, folderName) {
      if (idx === this.toggleFolder) {
        this.toggleFolder = null;
      } else {
        this.toggleFolder = idx;
        this.currentFolder = folderName;
        const folderFile = await findFolder({
          folder_name: folderName,
        });
        this.files_list = folderFile.data;
      }
    },
    async fileClick(idx, fileName) {
      this.apply = false;
      if (idx === this.toggleFile) {
        this.toggleFile = null;
      } else {
        this.toggleFile = idx;
        const file = await findFile({
          file_name: fileName,
        });
        this.selectFile = file.data;
        try {
          await checkConvert(this.selectFile.file_name);
          // this.apply = false;
        } catch (error) {
          console.error(error);
        }
      }
    },
    applyFile() {
      this.isLoading = true;
      try {
        // const file = await convertFile({
        //   file_name: this.selectFile.file_name,
        // });
        this.apply = true;
        const file_info = {
          file_name: this.selectFile.file_name,
          id: this.nodeId,
        };
        this.$store.commit("setWorkflowFile", file_info);
        this.$store.commit("shareWorkflowFile", this.nodeId);
      } catch (error) {
        console.error(error);
      }
      this.isLoading = false;
    },
  },
  async mounted() {
    const currentFile = this.$store.getters.getFileInfo(this.nodeId);
    await this.getFinder();
    if (currentFile.file !== "") {
      try {
        const file = await findFile({
          file_name: this.filterAndAddSuffix(currentFile.file),
        });
        this.selectFile = file.data;
        await checkConvert(this.filterAndAddSuffix(currentFile.file));
        this.apply = true;
      } catch (error) {
        console.error(error);
      }
    }
  },
  filters: {
    formatBytes(a, b) {
      if (a === 0) return "0 Bytes";
      const c = 1024;
      const d = b || 2;
      const e = ["Bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"];
      const f = Math.floor(Math.log(a) / Math.log(c));

      return parseFloat((a / Math.pow(c, f)).toFixed(d)) + " " + e[f];
    },
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
}

.fileUpload-form {
  width: 100%;
  height: 100%;
  position: relative;
  display: flex;
  align-items: flex-start;
}

.cloud-form {
  width: 45%;
  height: 95%;
  display: flex;
  flex-direction: column;
  /* align-items: stretch; */
  justify-content: space-evenly;
  margin: 1rem 0 1rem 1rem;
  padding: 8rem 3rem 3rem 3rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}

.folder__list {
  width: 45%;
  height: 95%;
  display: flex;
  flex-direction: column;
  margin: 1rem 0 1rem 1rem;
  padding: 1rem 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}

.folder__item {
  width: 100%;
  height: 5%;
  display: flex;
  cursor: pointer;
  color: #ffffff;
  margin-bottom: 0.5rem;
}

.file__item {
  width: 90%;
  height: 5%;
  margin-left: auto;
  display: flex;
  cursor: pointer;
  margin-bottom: 0.5rem;
  /* background: #ffffff; */
}

.folder__item:hover {
  background: rgb(231, 233, 238);
}

.toggleFolder {
  background: rgb(231, 233, 238);
}

.toggleFile {
  background: rgb(231, 233, 238);
}

.folder__item--col {
  width: 5rem;
  height: 100%;
  display: flex;
  align-items: center;
}

.folder__item--arrow {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
}

.folder__list--back {
  margin: 0 0 0.5rem 0.5rem;
  color: #2f2f2f;
  cursor: pointer;
  font-weight: 300;
}

.folder__list--back:hover {
  /* color: #000000; */
  font-weight: 400;
}

.folder__item--icon {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: contain;
  margin: 0 0.5rem;
}

.folder__name {
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: black;
}

.cloud-row {
  width: 100%;
  padding: 5% 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  text-decoration: none;
}

.cloud-row:first-child:hover {
  cursor: pointer;
}

.cloud-textbox {
  width: 17rem;
  height: 3rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  padding: 0 0 0 0rem;
}

.cloud-textbox__title {
  color: rgb(51, 51, 51);
  font-size: 2rem;
  font-weight: 400;
}

.cloud-textbox__desc {
  color: rgb(51, 51, 51);
  font-size: 1rem;
  font-weight: 200;
  padding: 0.2rem 0 0 0;
}

.cloud-textbox__directory {
  color: rgb(109, 158, 235);
  font-size: 1.2rem;
  font-weight: 400;
}

.cloud-textbox__directory__desc {
  color: rgb(51, 51, 51);
  font-size: 0.8rem;
  font-weight: 200;
  padding: 0.2rem 0 0 0;
}

.fileUpload {
  width: 55%;
  height: 100%;
  margin: 0rem 2rem 0rem 0;
  display: flex;
  flex-direction: column;
  align-items: center;
}

.form-row {
  width: 100%;
  height: 10%;
  position: relative;
}

.form-row:first-child {
  height: 75%;
  margin: 1rem;
  /* background: #3498db; */
  padding: 0;
}

.form-row:last-child {
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  height: 15%;
  margin: 0rem;
  padding: 0 2rem;
  box-sizing: border-box;
}

.form__name {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  padding: 0 2rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.4rem;
  line-height: 1.4rem;
  color: rgb(51, 51, 51);
}

.form__selectFile {
  width: 100%;
  height: 90%;
  display: flex;
  padding: 0 0 0 2rem;
}

.form__fileList {
  width: 84%;
  height: 100%;
  display: flex;
  padding: 2rem 1rem;
  flex-direction: column;
  align-items: center;
  overflow: hidden;
  border-radius: 1rem;
  background-color: rgb(255, 255, 255);
}

.fileList__item {
  width: 90%;
  padding: 0.5rem;
  background: #ffffff;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  margin-bottom: 0.5rem;
}

.fileList__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(53, 54, 58);
}

.fileList__text--blank {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(51, 51, 51, 0.5);
}

.form__button {
  width: 4rem;
  height: 4rem;
  display: flex;
  align-items: center;
  justify-content: center;
  margin: 1rem;
}

.form__addfile {
  width: 9rem;
  min-width: 9rem;
  height: 9rem;
  position: relative;
  cursor: pointer;
  background-color: rgb(80, 120, 251);
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 2rem;
  margin-top: -4rem;
}

.form__addfile:hover {
  background: rgb(96, 146, 246);
}

.form__button--plusicon {
  width: 1rem;
  height: 1rem;
  object-fit: cover;
  position: absolute;
  top: 2rem;
  right: 2rem;
}

.form__button--foldericon {
  width: 7rem;
  height: 7rem;
  object-fit: contain;
  position: absolute;
  top: 1rem;
  right: 1rem;
}

.form__input {
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
}

.form__info {
  width: 80%;
  height: 60%;
  display: flex;
  /* margin: auto; */
  top: 0;
  background: #ffffff;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
}

.form__info--type {
  width: 20%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}

.form__info--img {
  width: 2rem;
  height: 2rem;
  object-fit: cover;
}

.form__info--name,
.form__info--size {
  width: 80%;
  height: 100%;
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.9rem;
  line-height: 0.9rem;
  color: rgb(51, 51, 51);
  overflow-wrap: break-word;
  /* 줄바꿈 처리 */
  overflow: hidden;
  /* 넘치는 내용 숨기기 */
}

.form__info--blank {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(51, 51, 51, 0.5);
}

.form__button--apply {
  cursor: pointer;
  position: absolute;
  left: 78%;
  top: 7%;
  width: 20%;
  height: 46%;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  /* background: rgba(16, 83, 217, 0.377); */
  background: rgb(40, 84, 197);
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(244, 246, 251);
}

.form__button--apply:hover {
  background: rgb(75, 119, 209);
}

.activate {
  background: rgb(40, 197, 105);
}

.activate:hover {
  background: rgb(40, 197, 105);
  /* background: rgb(75, 232, 140); */
}

.loading-layout {
  position: absolute;
  top: 14px;
  right: -25px;
  display: flex;
  align-items: center;
  justify-content: center;
}

.loading-layout span {
  border: 4px solid #f3f3f3;
  border-top: 4px solid #3498db;
  border-radius: 50%;
  width: 20px;
  height: 20px;
  animation: spin 2s linear infinite;
}

/* @media (prefers-color-scheme: dark) {
  .form__button--foldericon {
    filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg)
      brightness(107%) contrast(101%);
  }
  .cloud-form {
    width: 45%;
    height: 95%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    margin: 1rem 0 1rem 1rem;
    padding: 3rem 0 3rem 1rem;
    border-radius: 1rem;
    box-sizing: border-box;
    background-color: rgb(41, 43, 48);
  }
  .folder__list {
    width: 45%;
    height: 95%;
    display: flex;
    flex-direction: column;
    margin: 1rem 0 1rem 1rem;
    padding: 3rem 1rem;
    border-radius: 1rem;
    box-sizing: border-box;
    background-color: rgb(41, 43, 48);
  }
  .folder__name {
    display: flex;
    align-items: center;
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 400;
    font-size: 1rem;
    line-height: 1rem;
    color: #ffffff;
  }
  .folder__item:hover {
    background: rgb(62, 64, 69);
  }
  .toggleFolder {
    background: rgb(62, 64, 69);
  }
  .toggleFile {
    background: rgb(62, 64, 69);
  }
  .folder__item--arrow {
    width: 1rem;
    height: 1rem;
    object-fit: contain;
    margin: 0 0.5rem;
    filter: invert(100%) sepia(0%) saturate(7487%) hue-rotate(164deg)
      brightness(106%) contrast(106%);
  }
  .form__addfile {
    width: 3rem;
    height: 3rem;
    position: relative;
    cursor: pointer;
    background-color: rgb(32, 34, 39);
    box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
    border-radius: 1rem;
  }
  .cloud-textbox__title {
    color: rgb(255, 255, 255);
    font-size: 1rem;
    font-weight: 400;
  }
  .cloud-textbox__desc {
    color: rgb(255, 255, 255);
    font-size: 1rem;
    font-weight: 200;
    padding: 0.2rem 0 0 0;
  }
  .form__name {
    width: 100%;
    height: 3rem;
    display: flex;
    align-items: center;
    padding: 0 2rem;
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 400;
    font-size: 1.4rem;
    line-height: 1.4rem;
    color: rgb(255, 255, 255);
  }
} */
</style>
