<template>
  <div class="layout" @click="ClickOut">
    <section class="folder">
      <div class="folder__header">
        <p class="folder__title left">Navigator</p>
        <div class="folder__createBtn right">
          <img class="folder__createBtn--icon" src="@/assets/add-folder.png" />
          <img class="folder__createBtn--icon" src="@/assets/refresh.png" />
        </div>
      </div>
      <ul class="folder__list">
        <li
          class="folder__item"
          v-for="(folder, idx) in folders_list"
          :key="idx"
          v-bind:class="{ toggleFolder: toggleFolder === idx }"
          @click="folderClick(idx, folder[0])"
        >
          <div class="folder__item--col">
            <img
              class="folder__item--icon"
              src="@/assets/arrow-bottom.png"
              v-if="toggleFolder === idx"
            />
            <img
              class="folder__item--icon"
              src="@/assets/arrow-right.png"
              v-else
            />
            <img
              class="folder__item--icon large"
              src="@/assets/open-folder.png"
              v-if="toggleFolder === idx"
            />
            <img
              class="folder__item--icon large"
              src="@/assets/folder.png"
              v-else
            />
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
        <input
          class="files__search"
          type="text"
          name="search"
          placeholder="search anything..."
        />
        <div class="header__column right">
          <label class="files__button">
            <img
              class="files__button--icon"
              src="@/assets/upload-file-black.png"
            />
            <input
              class="files__input"
              type="file"
              name="file"
              ref="selectFile"
              @change.prevent="uploadFile"
            />
          </label>
          <button class="files__button">
            <img class="files__button--icon" src="@/assets/delete.png" />
          </button>
          <button class="files__button right">
            <img class="files__button--icon" src="@/assets/setting.png" />
          </button>
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
          <tr
            class="files__item"
            v-for="(file, idx) in files_list"
            :key="idx"
            @contextmenu.prevent
            @click.right="RMouseClick($event, file.file_name, idx)"
            v-bind:class="{ select: R_Mouse_isActive }"
          >
            <td>{{ file.file_name | cutFromDotName }}</td>
            <td>{{ file.created_at | cutFromT }}</td>
            <td>{{ file.file_name | cutFromDotType }}</td>
            <td>{{ file.file_size | formatBytes }}</td>
          </tr>
        </tbody>
      </table>
      <ul
        ref="filesMenu"
        class="files_menu"
        v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }"
      >
        <li>view</li>
        <li>plot</li>
        <li>rename</li>
        <li @click="removeFile">delete</li>
      </ul>
    </main>
    <div class="message" v-bind:class="{ toggleMessage: !toggleMessage }">
      <p class="message__text">{{ messageContent }}</p>
      <p class="message__undo" @click="undoDeletion">undo</p>
      <img
        class="message__close"
        @click="toggleMessage = !toggleMessage"
        src="@/assets/close.png"
      />
    </div>
  </div>
</template>

<script>
import { uploadForm, getFiles, findFolder, deleteFile } from "@/api/index";

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
      toggleMessage: false,
      deletionTimer: null,
      messageContent: "",
      targetFile: null,
    };
  },

  methods: {
    RMouseClick(event, file_name, idx) {
      this.R_Mouse_isActive = false;
      this.xPosition = event.clientX + "px";
      this.yPosition = event.clientY - 55 + "px";
      this.R_Mouse_isActive = true;
      this.file_name = file_name;
      this.list_idx = idx;
      console.log(event);
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
        this.selectFile = new File(
          [this.$refs.selectFile.files[0]],
          `${this.currentFolder}_${this.$refs.selectFile.files[0].name}`
        );
        try {
          const form = new FormData();
          form.append("files", this.selectFile);
          const response = await uploadForm(form);
          console.log(response);
          const folderList = await findFolder({
            folder_name: this.currentFolder,
          });
          console.log(folderList.data);
          this.files_list = folderList.data;
        } catch (error) {
          console.error(error);
        }
      }
    },
    removeFile() {
      this.targetFile = this.files_list[this.list_idx];
      this.files_list.splice(this.list_idx, 1);
      console.log(this.targetFile);
      this.toggleMessage = true;
      // 10초 안에 toggleMessage가 false로 바뀌면 deleteFile 실행 안 함, 안 바뀌면 실행
      this.messageContent = `${this.targetFile.file_name} is deleted`;
      this.deletionTimer = setTimeout(async () => {
        try {
          const file = {
            file_name: this.file_name,
          };
          console.log(file);
          const targetFile = await deleteFile(file);
          console.log(targetFile);
        } catch (error) {
          console.error(error);
        }
        this.toggleMessage = false;
      }, 10000);
    },
    undoDeletion() {
      this.files_list.push(this.targetFile);
      this.toggleMessage = false;
      clearTimeout(this.deletionTimer);
    },
  },

  async mounted() {
    if (this.$route.query.doUploadFile) {
      this.$refs.selectFile.click();
    }
    try {
      const fileList = await getFiles();
      console.log(fileList.data);
      this.folders_list = fileList.data;
      const folderList = await findFolder({
        folder_name: this.currentFolder,
      });
      console.log(folderList.data);
      this.files_list = folderList.data;
    } catch (error) {
      console.error(error);
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
    cutFromT(value) {
      return value.split("T")[0];
    },
    cutFromDotName(value) {
      return value.split(".")[0];
    },
    cutFromDotType(value) {
      return value.split(".")[1];
    },
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  position: relative;
  overflow: hidden;
}
.folder {
  width: 25%;
  height: 100%;
  background: #ffffff;
  border-right: 1px solid #e1e1e1;
}
.folder__header {
  width: 90%;
  height: 10%;
  margin: auto;
  display: flex;
  align-items: center;
  position: relative;
}
.folder__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.25rem;
  line-height: 1.25rem;
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
}
.folder__list {
  width: 85%;
  height: 90%;
  margin: auto;
}
.folder__item {
  width: 100%;
  height: 5%;
  display: flex;
  cursor: pointer;
}
.folder__item:hover {
  background: rgb(204, 218, 245);
}
.toggleFolder {
  background: rgb(204, 218, 245);
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
  width: 33%;
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
  font-weight: 500;
  font-size: 1.25rem;
  line-height: 1.25rem;
}
.files__button {
  width: 2rem;
  height: 2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  background: #ffffff;
  margin-right: 1rem;
}
.files__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
}
.files__input {
  display: none;
}
.files__search {
  width: 16rem;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 0.5rem;
  outline-style: none;
  background: #f7f7f7;
}
.files__search:focus {
  box-shadow: rgba(0, 0, 0, 0.35) 0px 5px 15px;
  border: none;
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

.files_menu {
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
}
.files_menu.open {
  display: block;
  opacity: 1;
  position: absolute;
}
.files_menu > li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}
.files_menu > li:hover {
  background: #e5e5e5;
}

.message {
  width: 20rem;
  height: 3rem;
  display: flex;
  align-items: center;
  justify-content: space-around;
  position: absolute;
  bottom: 1rem;
  left: calc(50% - 10rem);
  background: rgba(0, 0, 0, 0.8);
  border-radius: 1rem;
  padding: 0 1rem;
}
.message__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: #ffffff;
}
.message__undo {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  color: #9196ff;
  text-decoration: underline;
  cursor: pointer;
}
.message__close {
  cursor: pointer;
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
  opacity: 0.5;
  filter: invert(100%) sepia(0%) saturate(0%) hue-rotate(0deg) brightness(100%)
    contrast(100%);
}
.toggleMessage {
  display: none;
}
</style>
