<template>
  <div class="layout">
    <div class="logo-component">
      <router-link class="logo" to="/main" v-if="!isWorkflowPage">
        <img class="logo__img" src="@/assets/cellcraft_logo.png" />
        <!-- <p class="logo__text">CELLCRAFT</p> -->
      </router-link>
      <router-link class="logo" to="/projects" v-else>
        <img class="logo__img" src="@/assets/cellcraft_logo.png" />
        <!-- <p class="logo__text">CELLCRAFT</p> -->
      </router-link>
    </div>
    <nav class="menu-component" v-if="!isWorkflowPage && !isUserLogin">
      <ul class="header-menu"></ul>
    </nav>
    <nav class="menu-component" v-else-if="!isWorkflowPage && isUserLogin">
      <ul class="header-menu">
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/projects">
            Workflows
          </router-link>
        </li>
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/files">
            Files
          </router-link>
        </li>
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/datasets">
            Datasets
          </router-link>
        </li>
        <!-- <li class="header-menu__item">
          <router-link class="header-menu__link" to="/plugins">
            Plugins
          </router-link>
        </li> -->
        <!-- <li class="header-menu__item" v-if="isSuperUser"> -->
        <li class="header-menu__item" v-if="isSuperUser">
          <router-link class="header-menu__link" to="/admin">
            Admin
          </router-link>
        </li>
      </ul>
    </nav>
    <div class="workflow-title" v-else>
      <p class="workflow-title__text" v-if="!activeInput" @click="editTitle">
        {{ setTitle }}
      </p>
      <input
        type="text"
        v-model="title"
        class="workflow-title__input"
        v-else
        @keydown.enter="applyTitle"
        @blur="applyTitle"
      />
    </div>
    <div class="login-component" v-if="isUserLogin">
      <p class="login__link" @click="logoutUser">Sign Out</p>
    </div>
    <div class="login-component" v-else>
      <router-link class="login__link" to="/login"> Sign In </router-link>
    </div>
  </div>
</template>

<script>
import { deleteCookie } from "@/utils/cookies";
import { getUser } from "@/api/index";

export default {
  data() {
    return {
      isWorkflowPage: false,
      title: this.$store.getters.getTitle,
      activeInput: false,
    };
  },
  computed: {
    isUserLogin() {
      return this.$store.getters.isLogin;
    },
    isSuperUser() {
      return this.$store.getters.isSuperUser;
    },
    setTitle() {
      return this.$store.getters.getTitle;
    },
  },
  methods: {
    redirect() {
      this.$router.push("/");
    },
    logoutUser() {
      this.$store.commit("clearToken");
      deleteCookie();
      this.$router.push("/");
    },
    editTitle() {
      this.activeInput = true;
    },
    applyTitle() {
      this.activeInput = false;
      this.$store.commit("setTitle", this.title);
    },
    async getUserInfo() {
      const response = await getUser();
      this.$store.commit("setUserInfo", response.data);
      console.log(response.data);
    },
  },
  watch: {
    $route(to) {
      this.isWorkflowPage = to.path.includes("/workflow");
      this.title = this.$store.getters.getTitle;
      this.getUserInfo();
    },
  },
};
</script>

<style scoped>
a {
  text-decoration: none;
  color: black;
}
.layout {
  height: 100%;
  width: 95rem;
  max-width: calc(100% - 3rem);
  margin-right: auto;
  margin-left: auto;
  display: flex;
  align-items: center;
  overflow: hidden;
}
.logo-component {
  width: 25%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: left;
  /* transform: translateY(3px); */
}
.logo {
  width: 2.25rem;
  height: 2.25rem;
  display: flex;
  align-items: center;
  justify-content: left;
}
.logo__img {
  width: 2.25rem;
  height: 2.25rem;
  object-fit: contain;
}
.logo__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 2rem;
  line-height: 2rem;
  text-decoration: none;
}
.menu-component {
  width: 50%;
  height: 100%;
}
.header-menu {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.header-menu__item {
  padding: 0 10%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  position: relative;
}
.header-menu__link {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(248, 250, 255);
  position: relative;
  transform: translateY(1px);
}
.header-menu__link:after {
  display: block;
  content: "";
  border-bottom: solid 2px rgb(248, 250, 255);
  transform: scaleX(0) translateY(4px);
  transition: transform 250ms ease-in-out;
}
.header-menu__link:hover:after {
  transform: scaleX(1) translateY(4px);
}
.header-menu__link.fromRight:after {
  transform-origin: 100% 50%;
}
.header-menu__link.fromLeft:after {
  transform-origin: 0% 50%;
}
.workflow-title {
  width: 50%;
  height: 100%;
  display: flex;
  justify-content: center;
  align-items: center;
}
.workflow-title__text {
  color: white;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.5rem;
  line-height: 1.5rem;
  text-decoration: none;
}
.workflow-title__input {
  color: white;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.5rem;
  line-height: 1.5rem;
  text-decoration: none;
  text-align: center;
  background: rgba(0, 0, 0, 0);
  border: none;
}
.workflow-title__input:focus {
  outline: none;
  box-shadow: none;
}
.login-component {
  width: 25%;
  height: 100%;
  display: flex;
  justify-content: right;
  align-items: center;
  /* transform: translateY(3px); */
}
.login__link {
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1.2rem;
  color: rgb(241, 243, 248);
}

/* @media (prefers-color-scheme: dark) {
  .header-menu__link {
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 400;
    font-size: 1rem;
    line-height: 1rem;
    color: rgb(248, 250, 255);
    position: relative;
    transform: translateY(4px);
  }
  .header-menu__link:after {
    display: block;
    content: "";
    border-bottom: solid 2px rgb(248, 250, 255);
    transform: scaleX(0) translateY(4px);
    transition: transform 250ms ease-in-out;
  }
  .login__link {
    display: flex;
    align-items: center;
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 600;
    font-size: 1.2rem;
    line-height: 1.2rem;
    color: rgb(241, 243, 248);
  }
} */
</style>
