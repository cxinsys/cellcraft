<template>
  <div class="layout">
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
        <li class="project-component">
          <div class="project__content">
            <img class="project__thumnail--icon" src="@/assets/create.png" />
          </div>
          <div class="project__info">
            <p class="project__title">New Workflow Project</p>
            <p class="project__date">Data Analysis Pipeline</p>
          </div>
        </li>
        <li class="project-component">
          <div class="project__content">
            <img
              class="project__thumnail"
              src="@/assets/workflow-template.png"
            />
          </div>
          <div class="project__info">
            <p class="project__title">Workflow Basics</p>
            <p class="project__date">Edited 2 hours ago</p>
          </div>
        </li>
        <li class="project-component">
          <div class="project__content">
            <img
              class="project__thumnail"
              src="@/assets/workflow-template.png"
            />
          </div>
          <div class="project__info">
            <p class="project__title">Workflow Basics2</p>
            <p class="project__date">Edited 2 hours ago</p>
          </div>
        </li>
      </ul>
    </div>
  </div>
</template>

<script>
import { getUser } from "@/api/index";

export default {
  data() {
    return {
      toggleMenu: true,
      profile: {},
    };
  },
  async mounted() {
    try {
      const profile = await getUser();
      console.log(profile.data);
      this.profile = profile.data;
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
  overflow: hidden;
  display: flex;
}
.side-panel {
  width: 20rem;
  height: 100%;
  border-right: 1px solid #e1e1e1;
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
}
.projects__recent,
.projects__template {
  width: 100%;
  height: 4rem;
  display: flex;
  align-items: center;
  cursor: pointer;
}
.projects__recent:hover,
.projects__template:hover {
  background: rgb(204, 218, 245);
}
.toggleMenu {
  background: rgb(204, 218, 245);
}
.projects__icon {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  margin-left: 1.5rem;
}
.projects__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.2rem;
  line-height: 1rem;
  padding-left: 1rem;
}
.project-view {
  width: calc(100% - 20rem);
  height: 100%;
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
  font-weight: 500;
  font-size: 1.4rem;
  line-height: 1rem;
  padding-left: 2rem;
}
.project-view__list {
  width: 100%;
  height: calc(100% - 5rem);
  display: flex;
}
.project-component {
  width: calc(33% - 3.85rem);
  height: 16rem;
  display: flex;
  flex-direction: column;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  margin: 3rem 0 0 3rem;
  cursor: pointer;
}
.project-component:hover {
  box-shadow: rgba(0, 0, 0, 0.1) 0px 4px 12px;
}
.project__content {
  width: 100%;
  height: 11rem;
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
</style>
