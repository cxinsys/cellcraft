<template>
  <div class="layout">
    <div class="first-line">
      <div class="first-line__left">
        <div class="header__text">
          Plugins
          <div class="header__desc">
            Plugin for Data Analysis Algorithm Extensions
          </div>
        </div>
      </div>
      <div class="first-line__right">
        <div class="add__button" @click="showPluginExtension = true">
          <img class="add__button--icon" src="@/assets/add_circle.png" />
          <h1>Add Plugin</h1>
        </div>
        <div class="search">
          <input
            type="text"
            v-model="searchTerm"
            placeholder="Search titles..."
          />
        </div>
      </div>
    </div>
    <PluginExtention v-if="showPluginExtension" @close="showPluginExtension = false" />
    <table>
      <tbody>
        <tr v-for="plugin in filteredPlugins" :key="plugin.id">
          <td>
            <div class="plugin-container">
              <div class="title-container">
                {{ plugin.title }}
              </div>
              <div class="description-container">
                {{ plugin.description }}
              </div>
              <div class="lastUpdated-container">
                Last Updated: {{ plugin.lastUpdated }}
              </div>
            </div>
            <div class="option-container">
              <div class="setting">
                <img
                  class="setting__button--icon"
                  src="@/assets/settings.png"
                />
              </div>
              <label class="switch">
                <input v-model="plugin.checked" type="checkbox" />
                <span class="slider round"></span>
              </label>
            </div>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>

<script>
import PluginExtention from "@/components/PluginExtention.vue";

export default {
  components: {
    PluginExtention,
  },
  data() {
    return {
      showPluginExtension: false,
      searchTerm: "",
      plugins: [
        {
          id: 1,
          title: "TENET",
          description:
            "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
          lastUpdated: "2024/03/18",
          checked: true,
        },
        {
          id: 2,
          title: "TENET TF",
          description:
            "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
          lastUpdated: "2024/03/18",
          checked: true,
        },
      ],
    };
  },
  computed: {
    filteredPlugins() {
      return this.plugins.filter((plugin) =>
        plugin.title.toLowerCase().includes(this.searchTerm.toLowerCase())
      );
    },
  },
  methods: {
    getCurrentDateString() {
      const today = new Date();
      const year = today.getFullYear();
      // 월은 0부터 시작하므로 1을 더해줍니다. 또한, 월과 일이 10보다 작을 때 앞에 '0'을 붙여줍니다.
      const month = String(today.getMonth() + 1).padStart(2, "0");
      const day = String(today.getDate()).padStart(2, "0");
      // YYYY/MM/DD 형식으로 문자열을 반환합니다.
      return `${year}/${month}/${day}`;
    },
  },
};
</script>

<style scoped>
.layout {
  padding: 10px 30px;
  overflow-y: auto;
}

table {
  width: 100%;
  height: 100%;
  border-collapse: separate;
  border-spacing: 5px;
  /* background-color: #c9c9c9; */
  transition: all 0.3s ease;
  border-radius: 15px;
  /* color: #ffffff; */
}

thead th,
td {
  padding: 10px;
  padding-left: 15px;
  text-align: left;
  border-radius: 10px;
  border: 1px solid #a8a8a8;
  /* box-shadow: 0px 4px 4px rgba(176, 169, 255, 0.25); */
}

th {
  text-transform: capitalize;
  background-color: #474747;
  color: #ffffff;
}

th:hover {
  background-color: #616161;
}

td {
  display: flex;
}

button {
  margin-right: 10px;
  color: black;
  padding: 5px;
  left: 10px;
  border-radius: 10px;
  background-color: #eaecff;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  text-transform: capitalize;
}

button:disabled {
  color: #ccc;
}

.sort-icon {
  color: rgb(199, 199, 199);
  font-weight: normal;
  font-size: small;
}

.first-line {
  width: calc(100% - 10px);
  margin: 1rem 5px;
  display: flex;
  align-items: center;
}

.first-line__left,
.first-line__right {
  width: calc(50% - 5px);
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.first-line__left {
  justify-content: left;
  align-items: end;
}

.header__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  /* padding-left: 2rem; */
  color: rgba(0, 0, 0, 0.8);
}

.header__desc {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(0, 0, 0, 0.5);
  display: flex;
  align-items: center;
  margin-top: 1rem;
}

.add__button {
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

.add__button:hover {
  cursor: pointer;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.add__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
}

.setting {
  width: 2rem;
  height: 2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
}

.setting__button--icon {
  width: 2rem;
  height: 2rem;
  object-fit: contain;
  opacity: 1;
}

.search {
  display: flex;
  align-items: center;
}

.search input {
  width: 300px;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 2rem;
  outline-style: none;
  background: #f7f7f7;
}

.search input:focus {
  border: 1px solid #bcbcbc;
}

#pageSize {
  padding: 2px;
  border-radius: 5px;
  border: 1px solid #ccc;
  margin-bottom: 5px;
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

.download-icon {
  margin: 0px 0px;
  width: 33px;
  height: 33px;
}

.plugin-container {
  width: calc(100% - 8rem);
}

.title-container {
  width: 100%;
  font-size: 1.4rem;
  align-items: center;
  display: flex;
  font-weight: 600;
  margin-top: 5px;
}

.description-container {
  width: 100%;
  font-size: 1.2rem;
  font-weight: 400;
  color: #474747;
  margin: 5px 2px;
  margin-bottom: 1rem;
}

.lastUpdated-container {
  width: 100%;
  font-size: 1rem;
  font-weight: 400;
  color: #474747;
  display: inline-block;
}

.option-container {
  width: 8rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
}

.switch {
  position: relative;
  display: inline-block;
  width: 60px;
  height: 34px;
  margin-top: 1rem;
}

.switch input {
  opacity: 0;
  width: 0;
  height: 0;
}

.slider {
  position: absolute;
  cursor: pointer;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background-color: #ccc;
  -webkit-transition: 0.4s;
  transition: 0.4s;
}

.slider:before {
  position: absolute;
  content: "";
  height: 26px;
  width: 26px;
  left: 4px;
  bottom: 4px;
  background-color: white;
  -webkit-transition: 0.4s;
  transition: 0.4s;
}

input:checked + .slider.w-color {
  background-color: #ccc;
}

input:checked + .slider.icon {
  background-color: #a37eba;
}

.slider.icon:before {
  background-color: #ffe05d;
}

.slider.icon:after {
  background-color: #e2df23;
}

input:checked + .slider {
  background-color: #2196f3;
}

input:checked + .slider:before {
  -webkit-transform: translateX(26px);
  -ms-transform: translateX(26px);
  transform: translateX(26px);
}

/* Rounded sliders */
.slider.round {
  border-radius: 34px;
}

.slider.round:before {
  border-radius: 50%;
}

.disabled-toggle {
  opacity: 0.5;
  cursor: default;
}
</style>
