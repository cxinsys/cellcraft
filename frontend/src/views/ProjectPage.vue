<template>
  <div class="layout" @click="ClickOut">
    <section class="side-panel">
      <div class="profile__header">
        <div class="profile__title">Projects</div>
      </div>
      <div class="profile">
        <img class="profile__img" src="@/assets/user.png" />
        <div class="profile__column">
          <p class="profile__username">{{ profile.username }}</p>
          <p class="profile__email">{{ profile.email }}</p>
        </div>
      </div>
      <div class="projects">
        <div class="projects__recent" v-bind:class="{ toggleMenu: toggleMenu }" @click="toggleMenu = true">
          <img class="projects__icon" src="@/assets/recent.png" />
          <p class="projects__text">Recents</p>
        </div>
      </div>
    </section>
    <div class="project-view">
      <header class="project-view__header">
        <p class="header__text">Recently Viewed</p>
      </header>
      <ul class="project-view__list">
        <li class="project-component" @click="createProject" @contextmenu.prevent>
          <div class="project__content">
            <img class="project__thumnail--icon" src="@/assets/create.png" />
          </div>
          <div class="project__info">
            <p class="project__title">New Workflow</p>
            <p class="project__date">Data Analysis Pipeline</p>
          </div>
        </li>
        <li class="project-component" v-for="(workflow, idx) in workflows" :key="idx" @contextmenu.prevent
          @click="openWorkflow(workflow.id)" @click.right="RMouseClick($event, workflow.id, idx)">
          <div class="project__content">
            <img class="project__thumnail" :src="workflow.thumbnail || require('@/assets/workflow-template2.png')
              " />
          </div>
          <div class="project__info">
            <p class="project__title">{{ workflow.title }}</p>
            <p class="project__date">{{ workflow.updated_at | updateTime }}</p>
          </div>
        </li>
      </ul>
      <ul ref="filesMenu" class="files_menu" v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }">
        <li @click="openWorkflow">Open</li>
        <li @click="confirmDelete">Delete</li>
      </ul>
    </div>
    <div class="message" v-bind:class="{ toggleMessage: !toggleMessage }">
      <img class="message__status" src="@/assets/succes.png" v-if="messageStatus === 'success'" />
      <img class="message__status" src="@/assets/error.png" v-else-if="messageStatus === 'error'" />
      <div class="message__box">
        <p class="message__text" v-for="(content, index) in filteredMessageContent" :key="index">
          {{ content }}
        </p>
      </div>
      <img class="message__close" @click="toggleMessage = !toggleMessage" src="@/assets/close.png" />
    </div>
    <!-- Modal for delete confirmation -->
    <div class="delete-modal" v-if="showDeleteModal">
      <div class="delete-modal__content">
        <p>Are you sure you want to delete this workflow ?</p>
        <ul class="delete-modal__buttons">
          <button @click="removeWorkflow">Yes</button>
          <button @click="showDeleteModal = !showDeleteModal">No</button>
        </ul>
        <img class="delete-modal__close" @click="showDeleteModal = !showDeleteModal" src="@/assets/close.png" />
      </div>
    </div>
    <div v-if="isSelectModalVisible" class="plugin-select-modal">
      <div class="plugin-select-modal__content">
        <!-- <span class="close" @click="closeSelectModal">&times;</span> -->
        <h2>Select a Plugin Template</h2>
        <ul>
          <li class="plugin-item" v-for="plugin in filteredPlugins" :key="plugin.id" @click="selectPlugin(plugin)">
            {{ plugin.name }}
            <span class="arrow">→</span>
          </li>
        </ul>
      </div>
    </div>
  </div>
</template>

<script>
import { getUser, getWorkflows, deleteWorkflow, getPlugins } from "@/api/index";

export default {
  data() {
    return {
      toggleMenu: true,
      profile: {},
      workflows: null,
      workflow_id: null,
      list_idx: null,
      R_Mouse_isActive: false,
      xPosition: 0,
      yPosition: 0,
      toggleMessage: false,
      targetWorkflow: null,
      plugins: [],
      isSelectModalVisible: false,
      showDeleteModal: false,
      messageStatus: "",
      messageContent: "",
    };
  },
  methods: {
    createProject() {
      this.isSelectModalVisible = true;
    },
    openWorkflow(workflow_id) {
      if (workflow_id) {
        this.$router.push({
          path: "/workflow",
          query: { workflow_id: workflow_id },
        });
      } else {
        this.$router.push({
          path: "/workflow",
          query: { workflow_id: this.workflow_id },
        });
      }
    },
    RMouseClick(event, id, idx) {
      this.R_Mouse_isActive = false;
      this.xPosition = event.clientX + "px";
      this.yPosition = event.clientY + "px";
      this.R_Mouse_isActive = true;
      this.workflow_id = id;
      this.list_idx = idx;
    },
    ClickOut() {
      this.R_Mouse_isActive = false;
    },
    confirmDelete() {
      this.targetWorkflow = this.workflows[this.list_idx];
      this.showDeleteModal = true;
    },
    async removeWorkflow() {
      try {
        const workflow = {
          id: this.workflow_id,
        };
        await deleteWorkflow(workflow);
        this.workflows.splice(this.list_idx, 1);
        this.showDeleteModal = false;
        this.setMessage("success", "Workflow deleted successfully");
      } catch (error) {
        console.error(error);
      }
    },
    closeSelectModal() {
      this.isSelectModalVisible = false;
    },
    selectPlugin(plugin) {
      const plugin_id = plugin.id;
      this.closeSelectModal();
      this.$router.push({
        path: "/workflow",
        query: { plugin_id: plugin_id },
      });
    },
    setMessage(status, content) {
      this.toggleMessage = true;
      this.messageStatus = status;
      this.messageContent = content;
      setTimeout(() => {
        this.toggleMessage = false;
      }, 5000);
    },
  },
  async mounted() {
    try {
      const [userResult, workflowsResult, pluginsResult] = await Promise.allSettled([
        getUser(),
        getWorkflows(),
        getPlugins()
      ]);

      // Handle getUser result
      if (userResult.status === "fulfilled") {
        const profile = userResult.value;
        this.profile = profile.data;
        const currentUser = profile.data.username;

        // Handle getWorkflows result
        if (workflowsResult.status === "fulfilled") {
          const workflows = workflowsResult.value;
          workflows.data.sort((a, b) => {
            return new Date(b.updated_at) - new Date(a.updated_at);
          });
          this.workflows = workflows.data;
        } else {
          console.error("Failed to fetch workflows:", workflowsResult.reason);
        }

        // Handle getPlugins result
        if (pluginsResult.status === "fulfilled") {
          const plugins = pluginsResult.value;
          console.log(plugins.data.plugins);
          this.plugins = plugins.data.plugins.map(plugin => {
            const userIncluded = plugin.users.some(user => user.username === currentUser);
            return {
              ...plugin,
              checked: userIncluded,
            };
          });
        } else {
          console.error("Failed to fetch plugins:", pluginsResult.reason);
        }
      } else {
        console.error("Failed to fetch user profile:", userResult.reason);
      }
    } catch (error) {
      console.error("Unexpected error:", error);
    }
  },
  filters: {
    updateTime(updated_at) {
      const time_diff = Date.now() - new Date(updated_at).getTime();
      const hours = Math.floor(time_diff / (1000 * 60 * 60));
      if (hours === 0) {
        return "Edited Recently";
      } else {
        return `Edited ${hours} hours ago`;
      }
    },
  },
  computed: {
    filteredPlugins() {
      return this.plugins.filter((plugin) => plugin.checked === true);
    },
    filteredMessageContent() {
      return this.messageContent.split(".").filter(Boolean);
    },
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  overflow: hidden;
  display: flex;
}

.side-panel {
  width: 18rem;
  height: 100%;
  border-right: 1px solid #aeaeae;
  background-color: rgb(201, 202, 203);
}

.profile {
  width: 100%;
  height: 4rem;
  display: flex;
  align-items: center;
  border-bottom: 1px solid #e1e1e1;
}

.profile__img {
  width: 2rem;
  height: 2rem;
  object-fit: cover;
  margin-left: 1.5rem;
  opacity: 0.8;
}

.profile__column {
  width: calc(100% - 3.5rem);
  height: 100%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  padding-left: 1rem;
}

.profile__username,
.profile__email {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.2rem;
  line-height: 1.2rem;
  color: rgba(0, 0, 0, 1);
}

.profile__email {
  font-weight: 400;
  font-size: 0.9rem;
  line-height: 0.9rem;
  color: rgba(0, 0, 0, 0.6);
}

.projects {
  width: 100%;
  display: flex;
  flex-direction: column;
  padding: 1rem 0;
}

.projects__recent,
.projects__template {
  width: calc(100% - 2rem);
  height: 3rem;
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  display: flex;
  align-items: center;
  cursor: pointer;
}

.projects__recent:hover,
.projects__template:hover {
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  background: rgb(176, 177, 178);
}

.toggleMenu {
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  background: rgb(176, 177, 178);
}

.projects__icon {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  margin-left: 1.5rem;
  color: rgba(0, 0, 0, 0.7);
}

.projects__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  padding-left: 1rem;
  color: rgba(0, 0, 0, 0.7);
}

.project-view {
  width: calc(100% - 20rem);
  height: 100%;
  overflow-y: auto;
  background-color: #ffffff;
}

.project-view__header {
  width: 100%;
  height: 5rem;
  display: flex;
  align-items: center;
  border-bottom: 1px solid #e1e1e1;
}

.header__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  padding-left: 2rem;
  color: rgba(0, 0, 0, 0.8);
}

.project-view__list {
  display: grid;
  grid-template-columns: repeat(auto-fill,
      minmax(230px, 1fr));
  /* 컴포넌트 최소 너비를 250px, 최대를 1fr로 설정 */
  padding: 3rem;
  gap: 3rem;
  /* for space around items */
}

.project-component {
  /* flex: 1 0 calc(33% - 6rem);
  height: 16rem;
  max-width: 19.5rem; */
  display: flex;
  flex-direction: column;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  cursor: pointer;
  box-shadow: rgba(0, 0, 0, 0.1) 0px 2px 2px;
}

.project-component:hover {
  box-shadow: rgba(0, 0, 0, 0.1) 0px 4px 6px;
}

.project__content {
  width: 100%;
  height: 10rem;
  display: flex;
  align-items: center;
  justify-content: center;
}

.project__thumnail {
  width: 100%;
  height: 100%;
  object-fit: cover;
  border-radius: 1rem 1rem 0 0;
}

.project__thumnail--icon {
  width: 4rem;
  height: 4rem;
  object-fit: cover;
}

.project__info {
  width: 100%;
  height: 5rem;
  border-top: 1px solid #e1e1e1;
  display: flex;
  flex-direction: column;
  justify-content: center;
}

.project__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.2rem;
  line-height: 1.3rem;
  padding: 0 0 0.5rem 1rem;
  color: rgba(0, 0, 0, 0.8);
}

.project__date {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.9rem;
  line-height: 1rem;
  padding-left: 1rem;
  opacity: 0.6;
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
  text-transform: capitalize;
}

.files_menu.open {
  display: block;
  opacity: 1;
  position: absolute;
}

.files_menu>li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}

.files_menu>li:hover {
  background: #e5e5e5;
}

.delete-modal {
  display: block;
  position: fixed;
  left: 0;
  top: 0;
  z-index: 1;
  height: 100%;
  width: 100%;
  overflow: auto;
  background-color: rgb(0, 0, 0);
  background-color: rgba(0, 0, 0, 0.4);
}

/* .delete-modal__content {
  background-color: #fefefe;
  margin: 15% auto;
  padding: 20px;
  border: 1px solid #888;
  width: 80%;
} */

.delete-modal__content {
  position: relative;
  width: 90%;
  max-width: 400px;
  margin: 25% auto;
  background: #FFF;
  border-radius: .25em .25em .4em .4em;
  text-align: center;
  box-shadow: 0 0 20px rgba(0, 0, 0, 0.2);
  -webkit-transform: translateY(-40px);
  -moz-transform: translateY(-40px);
  -ms-transform: translateY(-40px);
  -o-transform: translateY(-40px);
  transform: translateY(-40px);
  /* Force Hardware Acceleration in WebKit */
  -webkit-backface-visibility: hidden;
  -webkit-transition-property: -webkit-transform;
  -moz-transition-property: -moz-transform;
  transition-property: transform;
  -webkit-transition-duration: 0.3s;
  -moz-transition-duration: 0.3s;
  transition-duration: 0.3s;
}

.delete-modal__content p {
  padding: 2.5rem 1rem;
}

.delete-modal__content .delete-modal__buttons:after {
  content: "";
  display: table;
  clear: both;
}

.delete-modal__content .delete-modal__buttons button {
  float: left;
  width: 50%;
  border: 0;
  cursor: pointer;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1rem;
  line-height: 1.25rem;
  text-align: center;
  text-transform: uppercase;
  padding: 1rem 0.5rem;
}

.delete-modal__content .delete-modal__buttons a {
  display: block;
  height: 60px;
  line-height: 60px;
  text-transform: uppercase;
  color: #FFF;
  -webkit-transition: background-color 0.2s;
  -moz-transition: background-color 0.2s;
  transition: background-color 0.2s;
}

.delete-modal__content .delete-modal__buttons button:first-child {
  background: #fc7169;
  border-radius: 0 0 0 .25em;
}

.no-touch .delete-modal__content .delete-modal__buttons button:first-child:hover {
  background-color: #fc8982;
}

.delete-modal__content .delete-modal__buttons button:last-child {
  background: #b6bece;
  border-radius: 0 0 .25em 0;
}

.no-touch .delete-modal__content .delete-modal__buttons button:last-child:hover {
  background-color: #c5ccd8;
}

.delete-modal__content .delete-modal__close {
  position: absolute;
  top: 8px;
  right: 8px;
  width: 1rem;
  height: 1rem;
  cursor: pointer;
  object-fit: contain;
}

.delete-modal__content .delete-modal__close::before,
.delete-modal__content .delete-modal__close::after {
  content: '';
  position: absolute;
  top: 12px;
  width: 14px;
  height: 3px;
  background-color: #8f9cb5;
}

.delete-modal__content .delete-modal__close::before {
  -webkit-transform: rotate(45deg);
  -moz-transform: rotate(45deg);
  -ms-transform: rotate(45deg);
  -o-transform: rotate(45deg);
  transform: rotate(45deg);
  left: 8px;
}

.delete-modal__content .delete-modal__close::after {
  -webkit-transform: rotate(-45deg);
  -moz-transform: rotate(-45deg);
  -ms-transform: rotate(-45deg);
  -o-transform: rotate(-45deg);
  transform: rotate(-45deg);
  right: 8px;
}

.message {
  width: 30rem;
  height: 6rem;
  display: flex;
  align-items: center;
  justify-content: space-around;
  position: absolute;
  bottom: 5rem;
  left: calc(50% - 16rem);
  background: rgba(0, 0, 0, 0.8);
  border-radius: 1rem;
  padding: 0 1rem;
}

.message__status {
  width: 2rem;
  height: 2rem;
  object-fit: contain;
  margin: 0 1rem;
}

.message__box {
  display: flex;
  flex-direction: column;
  justify-content: center;
  width: 100%;
  height: 100%;
}

.message__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.1rem;
  line-height: 1.4rem;
  color: #ffffff;
}

.message__close {
  cursor: pointer;
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
  opacity: 0.5;
  filter: invert(100%) sepia(0%) saturate(0%) hue-rotate(0deg) brightness(100%) contrast(100%);
}

.toggleMessage {
  display: none;
}

.profile__header {
  width: 90%;
  height: 10%;
  margin: auto;
  display: flex;
  align-items: center;
  position: relative;
  color: rgba(0, 0, 0, 0.8);
}

.profile__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.7rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
}

.plugin-select-modal {
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(0, 0, 0, 0.5);
  display: flex;
  justify-content: center;
  align-items: center;
}

.plugin-select-modal__content {
  background-color: white;
  border-radius: 10px;
  padding: 1.5rem 2rem;
  position: relative;
}

.plugin-select-modal__content h2 {
  width: 15rem;
  margin-bottom: 1rem;
}

.plugin-select-modal__content ul {
  list-style: none;
  padding: 0;
}

.plugin-select-modal__content li {
  padding: 10px;
  cursor: pointer;
}

.plugin-select-modal__content li:hover {
  background-color: #f5f5f5;
}

.plugin-item {
  position: relative;
  padding-right: 1rem;
}

.plugin-item .arrow {
  font-size: 1.1rem;
  position: absolute;
  right: 10px;
  opacity: 0;
  transition: all 0.5s;
}

.plugin-item:hover .arrow {
  opacity: 1;
  right: 5px;
}

/* .close {
  position: absolute;
  top: 1rem;
  right: 1rem;
  color: rgb(255, 90, 90);
  font-size: 28px;
  font-weight: bold;
}

.close:hover,
.close:focus {
  color: red;
  text-decoration: none;
  cursor: pointer;
} */
</style>
