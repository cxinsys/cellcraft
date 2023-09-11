<template>
  <div class="layout">
    <div class="side_nav" :style="{ top: sideNavTop }">
      <div
        v-for="button in buttons"
        :key="button.id"
        :class="['nav_button', { active: isActiveButton(button) }]"
        @click="navigateTo(button.id, button.index)"
      >
        <img class="button_img" :src="button.img" draggable="false" />
      </div>
    </div>
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
          label: "버튼 1",
          img: require("@/assets/userAdmin.png"),
        },
        {
          id: "admin/job",
          index: 1,
          label: "버튼 2",
          img: require("@/assets/control_jobs.png"),
        },
        {
          id: "admin/dataset",
          index: 2,
          label: "버튼 3",
          img: require("@/assets/datasetAdmin.png"),
        },
        {
          id: "admin/algorithm",
          index: 3,
          label: "버튼 4",
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
  margin-top: 60px;
  min-height: 100vh;
  height: 100%;
  display: flex;
}
.side_nav {
  left: 0.3rem;
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
  justify-content: space-evenly;
}
.nav_button {
  display: flex;
  width: 3rem;
  height: 3rem;
  padding: 0.1rem;
  margin: 0.1rem;
  margin-bottom: 0.5rem;
  border-radius: 20%;
  background-color: rgba(255, 255, 255, 0.1);
  transition: background-color 0.3s;
  justify-content: center;
  align-items: center;
}
.nav_button:hover {
  background-color: lightgray;
  cursor: pointer;
}

.nav_button.active {
  background-color: rgba(211, 211, 211, 0.612);
  transition: background-color 0.3s;
}
.nav_button.active:hover {
  background-color: lightgray;
  cursor: pointer;
}
.admin-component {
  margin: 0rem 4.5rem;
  width: calc(100vw - 3.5rem);
}
.admin_view {
  overflow-y: auto;
}

.button_img {
  width: 2rem;
  height: 2em;
}
</style>
