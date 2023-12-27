<template>
  <div class="layout">
    <section class="side_nav" :style="{ top: sideNavTop }">
      <div class="admin__header">
        <div class="admin__title">Admin</div>
      </div>
      <div
        v-for="button in buttons"
        :key="button.id"
        :class="['nav_button', { active: isActiveButton(button) }]"
        @click="navigateTo(button.id, button.index)"
      >
        <img class="button_img" :src="button.img" draggable="false" />
        <p class="button_text">{{ button.label }}</p>
      </div>
    </section>
    <div class="admin-component">
      <router-view></router-view>
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      buttons: [
        {
          id: "admin/user",
          index: 0,
          label: "Users",
          img: require("@/assets/userAdmin.png"),
        },
        {
          id: "admin/job",
          index: 1,
          label: "Jobs",
          img: require("@/assets/control_jobs.png"),
        },
        {
          id: "admin/dataset",
          index: 2,
          label: "Datasets",
          img: require("@/assets/datasetAdmin.png"),
        },
        {
          id: "admin/algorithm",
          index: 3,
          label: "Algorithms",
          img: require("@/assets/algorithm3.png"),
        },
      ],
      lastIndex: 0,
      sideNavTop: "20%",
    };
  },
  mounted() {
    window.addEventListener("scroll", this.handleScroll);
  },
  beforeDestroy() {
    window.removeEventListener("scroll", this.handleScroll);
  },
  methods: {
    navigateTo(id, index) {
      this.$router.push(`/${id}`);
      this.lastIndex = index;
    },
    handleScroll() {
      const navOffset = document.querySelector(".layout").offsetTop;
      const scrollTop =
        window.pageYOffset || document.documentElement.scrollTop;
      this.sideNavTop = navOffset <= scrollTop ? "20%" : "20%";
    },
  },
  computed: {
    isActiveButton() {
      return (button) => button.index === this.lastIndex;
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
.side_nav {
  /* left: 0.3rem;
  position: fixed;
  display: flex;
  width: 3.6rem;
  top: 20%;
  border-radius: 0.5rem;
  height: 16rem;
  padding-top: 10px;
  background-color: rgb(76, 76, 76);
  flex-direction: column;
  align-items: center;
  justify-content: space-evenly; */
  width: 18rem;
  height: 100%;
  border-right: 1px solid #aeaeae;
  background-color: rgb(201, 202, 203);
}
.nav_button {
  width: clac(100% - 2rem);
  height: 3rem;
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  display: flex;
  align-items: center;
  cursor: pointer;
  /* display: flex;
  width: 3rem;
  height: 3rem;
  padding: 0.1rem;
  margin: 0.1rem;
  margin-bottom: 0.5rem;
  border-radius: 20%;
  background-color: rgba(255, 255, 255, 0.1);
  transition: background-color 0.3s;
  justify-content: center;
  align-items: center; */
}
.nav_button:hover {
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  background: rgb(176, 177, 178);
  /* background-color: lightgray;
  cursor: pointer; */
}

.nav_button.active {
  margin: 0.1rem 1rem;
  border-radius: 1rem;
  background: rgb(176, 177, 178);
  /* background-color: rgba(211, 211, 211, 0.612);
  transition: background-color 0.3s; */
}
/* .nav_button.active:hover {
  background-color: lightgray;
  cursor: pointer;
} */
.admin-component {
  margin: 1rem 0rem;
  width: calc(100% - 18rem);
}
.admin_view {
  overflow-y: auto;
}

.button_img {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: contain;
  margin-left: 1.5rem;
  color: rgba(0, 0, 0, 0.7);
  filter: invert(100%);
}

.admin__header {
  width: 90%;
  height: 10%;
  margin: auto;
  display: flex;
  align-items: center;
  position: relative;
  color: rgba(0, 0, 0, 0.8);
}
.admin__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.7rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
}

.button_text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  padding-left: 1rem;
  color: rgba(0, 0, 0, 0.7);
}
</style>
