<template>
  <div class="layout" @click="ClickOut">
    <section class="side-panel">
      <div class="profile">
        <img class="profile__img" src="@/assets/user.png" />
        <div class="profile__column">
          <p class="profile__username">{{ profile.username }}</p>
          <p class="profile__email">{{ profile.email }}</p>
        </div>
      </div>
      <div class="projects">
        <div
          class="projects__recent"
          v-bind:class="{ toggleMenu: toggleMenu }"
          @click="toggleMenu = true"
        >
          <img class="projects__icon" src="@/assets/recent.png" />
          <p class="projects__text">Recents</p>
        </div>
        <div
          class="projects__template"
          v-bind:class="{ toggleMenu: !toggleMenu }"
          @click="toggleMenu = false"
        >
          <img class="projects__icon" src="@/assets/template.png" />
          <p class="projects__text">Templates</p>
        </div>
      </div>
    </section>
    <div class="project-view">
      <header class="project-view__header">
        <p class="header__text">Recently Viewed</p>
      </header>
      <ul class="project-view__list">
        <li
          class="project-component"
          @click="createProject"
          @contextmenu.prevent
        >
          <div class="project__content">
            <img class="project__thumnail--icon" src="@/assets/create.png" />
          </div>
          <div class="project__info">
            <p class="project__title">New Workflow Project</p>
            <p class="project__date">Data Analysis Pipeline</p>
          </div>
        </li>
        <!-- <li class="project-component" @contextmenu.prevent>
          <div class="project__content">
            <img
              class="project__thumnail"
              src="@/assets/workflow-template2.png"
            />
          </div>
          <div class="project__info">
            <p class="project__title">Workflow Basics</p>
            <p class="project__date">Just Template</p>
          </div>
        </li>
        <li class="project-component" @contextmenu.prevent>
          <div class="project__content">
            <img
              class="project__thumnail"
              src="@/assets/workflow-template2.png"
            />
          </div>
          <div class="project__info">
            <p class="project__title">Workflow Basics2</p>
            <p class="project__date">Just Template</p>
          </div>
        </li> -->
        <li
          class="project-component"
          v-for="(workflow, idx) in workflows"
          :key="idx"
          @contextmenu.prevent
          @click="openWorkflow(workflow.id)"
          @click.right="RMouseClick($event, workflow.id, idx)"
        >
          <div class="project__content">
            <!-- <img
              class="project__thumnail"
              src="@/assets/workflow-template2.png"
            /> -->
            <img
              class="project__thumnail"
              :src="
                workflow.thumbnail || require('@/assets/workflow-template2.png')
              "
            />
          </div>
          <div class="project__info">
            <p class="project__title">{{ workflow.title }}</p>
            <p class="project__date">{{ workflow.updated_at | updateTime }}</p>
          </div>
        </li>
      </ul>
      <ul
        ref="filesMenu"
        class="files_menu"
        v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }"
      >
        <li @click="openWorkflow">Open</li>
        <li>Copy Link</li>
        <li>Rename</li>
        <li @click="removeWorkflow">Delete</li>
      </ul>
    </div>
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
import { getUser, getWorkflows, deleteWorkflow } from "@/api/index";

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
      deletionTimer: null,
      messageContent: "",
      targetWorkflow: null,
    };
  },
  methods: {
    createProject() {
      this.$router.push("/workflow");
    },
    openWorkflow(id) {
      if (id) {
        this.$router.push({
          path: "/workflow",
          query: { id: id },
        });
      } else {
        this.$router.push({
          path: "/workflow",
          query: { id: this.workflow_id },
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
      console.log(event);
    },
    ClickOut() {
      this.R_Mouse_isActive = false;
    },
    async removeWorkflow() {
      this.targetWorkflow = this.workflows[this.list_idx];
      this.workflows.splice(this.list_idx, 1);
      console.log(this.targetWorkflow);
      this.toggleMessage = true;
      // 10초 안에 toggleMessage가 false로 바뀌면 deleteFile 실행 안 함, 안 바뀌면 실행
      this.messageContent = `${this.targetWorkflow.title} is deleted`;
      this.deletionTimer = setTimeout(async () => {
        try {
          const workflow = {
            id: this.workflow_id,
          };
          const targetWorkflow = await deleteWorkflow(workflow);
          console.log(targetWorkflow);
          this.workflows.splice(this.list_idx, 1);
        } catch (error) {
          console.error(error);
        }
        this.toggleMessage = false;
      }, 3000);
    },
    undoDeletion() {
      this.workflows.push(this.targetWorkflow);
      this.toggleMessage = false;
      clearTimeout(this.deletionTimer);
    },
  },
  async mounted() {
    try {
      const profile = await getUser();
      console.log(profile.data);
      this.profile = profile.data;
      const workflows = await getWorkflows();
      console.log(1);
      console.log(workflows.data);
      this.workflows = workflows.data;
    } catch (error) {
      console.error(error);
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
  width: 14rem;
  height: 100%;
  border-right: 1px solid #aeaeae;
  background-color: rgb(201, 202, 203);
}
.profile {
  width: 100%;
  height: 8rem;
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
  font-size: 1rem;
  line-height: 1rem;
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
  width: clac(100% - 2rem);
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
  grid-template-columns: repeat(
    auto-fill,
    minmax(230px, 1fr)
  ); /* 컴포넌트 최소 너비를 250px, 최대를 1fr로 설정 */
  padding: 3rem;
  gap: 3rem; /* for space around items */
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
