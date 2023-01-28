<template>
  <div class="layout">
    <div class="logo-component">
      <router-link class="logo" to="/main">
        <img class="logo__img" src="@/assets/cellcraft_logo.png" />
        <!-- <p class="logo__text">CELLCRAFT</p> -->
      </router-link>
    </div>
    <nav class="menu-component">
      <ul class="header-menu">
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/workflow">
            Workflows
          </router-link>
        </li>
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/files">
            Files
          </router-link>
        </li>
        <li class="header-menu__item">
          <router-link class="header-menu__link" to="/main">
            Datasets
          </router-link>
        </li>
      </ul>
    </nav>
    <div class="login-component">
      <router-link class="login__link" to="/login"> LOGIN </router-link>
    </div>
  </div>
</template>

<script>
import { getUser } from "@/api/index";
import { deleteCookie } from "@/utils/cookies";
// import Profile from "@/components/profile";

export default {
  // components: { Profile },
  data() {
    return {
      profile: {
        username: null,
        email: null,
        password: "********",
      },
    };
  },
  computed: {
    isUserLogin() {
      console.log(this.$store.getters.isLogin);
      return this.$store.getters.isLogin;
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
    async getProfile() {
      this.modal = true;
      const userInfo = await getUser();
      this.$store.commit("setUserInfo", userInfo.data);
      console.log(userInfo.data);
      this.profile.username = userInfo.data.username;
      this.profile.email = userInfo.data.email;
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
}
.logo-component {
  width: 25%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: left;
}
.logo {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: left;
}
.logo__img {
  width: 2.5rem;
  height: 2.5rem;
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
  color: rgb(70, 70, 70);
  position: relative;
  transform: translateY(1px);
}
.header-menu__link:after {
  display: block;
  content: "";
  border-bottom: solid 3px #575757;
  transform: scaleX(0) translateY(10px);
  transition: transform 250ms ease-in-out;
}
.header-menu__link:hover:after {
  transform: scaleX(1) translateY(10px);
}
.header-menu__link.fromRight:after {
  transform-origin: 100% 50%;
}
.header-menu__link.fromLeft:after {
  transform-origin: 0% 50%;
}

.login-component {
  width: 25%;
  height: 100%;
  display: flex;
  justify-content: right;
  align-items: center;
}
.login__link {
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.2rem;
  line-height: 1.2rem;
  color: rgba(51, 51, 51, 0.6);
}
</style>
