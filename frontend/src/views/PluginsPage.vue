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
        <div class="add__button" @click="addPluginExtension">
          <img class="add__button--icon" src="@/assets/add_circle.png" />
          <h1>Add Plugin</h1>
        </div>
        <!-- <div class="build-all__button" @click="buildAllPlugins">
          <img class="build-all__button--icon" src="@/assets/add_circle.png" />
          <h1>Build All</h1>
        </div> -->
        <div class="search">
          <input type="text" v-model="searchTerm" placeholder="Search titles..." />
        </div>
      </div>
    </div>
    <PluginExtention v-if="showPluginExtension" @close="closePluginExtension" :editName="selectedPlugin.name"
      :editDescription="selectedPlugin.description" :editDependencies="selectedPlugin.dependencies"
      :editDrawflow="selectedPlugin.drawflow" :editRules="selectedPlugin.rules" />
    <table>
      <tbody>
        <tr v-for="plugin in filteredPlugins" :key="plugin.id">
          <td>
            <div class="plugin-container">
              <div class="title-container">
                {{ plugin.name }}
              </div>
              <div class="description-container">
                {{ plugin.description }}
              </div>
              <div class="lastUpdated-container">
                Last Updated: {{ plugin.updated_at.split('T')[0] }}
              </div>
            </div>
            <div class="option-container">
              <div class="setting" @click="editPluginExtension(plugin)">
                <img class="setting__button--icon" src="@/assets/settings.png" />
              </div>
              <button class="build-button" @click="handleBuildPlugin(plugin)"
                :disabled="plugin.building || plugin.imageExists"
                :class="{ 'building': plugin.building, 'image-exists': plugin.imageExists }">
                <span v-if="plugin.building">Building...</span>
                <span v-else-if="plugin.imageExists">Built</span>
                <span v-else>Build</span>
              </button>
              <label class="switch">
                <input v-model="plugin.checked" type="checkbox" @change="handlePluginAssociate(plugin)"
                  :disabled="isCheckboxDisabled" />
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
import { getUser, getPlugins, associatePlugin, dissociatePlugin, buildPluginDocker, checkPluginImage } from "@/api/index";
import PluginExtention from "@/components/PluginExtention.vue";

export default {
  components: {
    PluginExtention,
  },
  data() {
    return {
      showPluginExtension: false,
      isCheckboxDisabled: false,
      searchTerm: "",
      plugins: [
        // {
        //   id: 1,
        //   name: "TENET",
        //   description:
        //     "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
        //   lastUpdated: "2024/03/18",
        //   checked: true,
        // },
        // {
        //   id: 2,
        //   name: "TENET TF",
        //   description:
        //     "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
        //   lastUpdated: "2024/03/18",
        //   checked: true,
        // },
      ],
      profile: {},
      selectedPlugin: {
        name: "",
        description: "",
        dependencies: {},
        drawflow: {},
        rules: [],
      },
    };
  },
  async mounted() {
    try {
      await this.getUserAssociatePlugins();
    } catch (error) {
      console.error(error);
    }
  },
  computed: {
    filteredPlugins() {
      return this.plugins.filter((plugin) =>
        plugin.name.toLowerCase().includes(this.searchTerm.toLowerCase())
      );
    },
  },
  methods: {
    async closePluginExtension() {
      this.showPluginExtension = false;
      // extension 완료했으니, 다시 plugin list를 불러옵니다.
      try {
        await this.getUserAssociatePlugins();
      } catch (error) {
        console.error(error);
      }
    },
    getCurrentDateString() {
      const today = new Date();
      const year = today.getFullYear();
      // 월은 0부터 시작하므로 1을 더해줍니다. 또한, 월과 일이 10보다 작을 때 앞에 '0'을 붙여줍니다.
      const month = String(today.getMonth() + 1).padStart(2, "0");
      const day = String(today.getDate()).padStart(2, "0");
      // YYYY/MM/DD 형식으로 문자열을 반환합니다.
      return `${year}/${month}/${day}`;
    },
    addPluginExtension() {
      this.selectedPlugin = {
        name: "",
        description: "",
        dependencies: {},
        drawflow: {},
        rules: [],
      };
      this.showPluginExtension = true;
    },
    editPluginExtension(plugin) {
      this.selectedPlugin = plugin;
      this.showPluginExtension = true;
    },
    async getUserAssociatePlugins() {
      try {
        const profile = await getUser();
        this.profile = profile.data;

        const plugins = await getPlugins();
        console.log(plugins.data.plugins);
        const currentUser = this.profile.username;

        this.plugins = plugins.data.plugins.map(plugin => {
          const userIncluded = plugin.users.some(user => user.username === currentUser);
          return {
            ...plugin,
            checked: userIncluded,
            building: false,
            imageExists: false,
          };
        });

        // 각 플러그인의 이미지 존재 여부 확인
        await this.checkAllPluginImages();
      } catch (error) {
        console.error(error);
      }
    },
    async handlePluginAssociate(plugin) {
      const pluginId = parseInt(plugin.id);

      try {
        let result;
        if (plugin.checked) {
          result = await associatePlugin(pluginId);
        } else {
          result = await dissociatePlugin(pluginId);
        }
        console.log(result.data);

        // 1초 동안 체크박스 비활성화
        this.isCheckboxDisabled = true;
        setTimeout(() => {
          this.isCheckboxDisabled = false;
        }, 1000);
      } catch (error) {
        console.error('Error associating/disassociating plugin:', error);
      }
    },
    async checkAllPluginImages() {
      for (let plugin of this.plugins) {
        try {
          const result = await checkPluginImage(plugin.name);
          plugin.imageExists = result.data.image_exists;
        } catch (error) {
          console.error(`Error checking image for plugin ${plugin.name}:`, error);
          plugin.imageExists = false;
        }
      }
    },
    async handleBuildPlugin(plugin) {
      if (plugin.building || plugin.imageExists) {
        return;
      }

      try {
        plugin.building = true;
        const result = await buildPluginDocker(plugin.name, false); // 기존 플러그인은 기본적으로 GPU 비활성화

        // 빌드 성공 후 이미지 존재 여부 다시 확인
        setTimeout(async () => {
          try {
            const checkResult = await checkPluginImage(plugin.name);
            plugin.imageExists = checkResult.data.image_exists;
          } catch (error) {
            console.error(`Error checking image after build for plugin ${plugin.name}:`, error);
          }
          plugin.building = false;
        }, 1000);

        console.log('Build result:', result.data);
        alert(`Plugin ${plugin.name} built successfully!`);
      } catch (error) {
        console.error(`Error building plugin ${plugin.name}:`, error);
        plugin.building = false;

        // 에러 메시지 표시
        let errorMessage = `Failed to build plugin ${plugin.name}`;
        if (error.response && error.response.data && error.response.data.detail) {
          if (typeof error.response.data.detail === 'string') {
            errorMessage += `: ${error.response.data.detail}`;
          } else if (error.response.data.detail.message) {
            errorMessage += `: ${error.response.data.detail.message}`;
          }
        }
        alert(errorMessage);
      }
    },
    async buildAllPlugins() {
      try {
        // 빌드가 필요한 플러그인들만 필터링 (building이 아니고 imageExists가 false인 것들)
        const pluginsToBuild = this.plugins.filter(plugin => !plugin.building && !plugin.imageExists);

        if (pluginsToBuild.length === 0) {
          alert('빌드할 플러그인이 없습니다. 모든 플러그인이 이미 빌드되었거나 빌드 중입니다.');
          return;
        }

        if (!confirm(`${pluginsToBuild.length}개의 플러그인을 빌드하시겠습니까?`)) {
          return;
        }

        const buildPromises = pluginsToBuild.map(plugin => {
          plugin.building = true;
          return buildPluginDocker(plugin.name, false) // 기존 플러그인은 기본적으로 GPU 비활성화
            .then(result => {
              console.log(`Build result for ${plugin.name}:`, result.data);
              setTimeout(async () => {
                try {
                  const checkResult = await checkPluginImage(plugin.name);
                  plugin.imageExists = checkResult.data.image_exists;
                } catch (error) {
                  console.error(`Error checking image after build for plugin ${plugin.name}:`, error);
                }
                plugin.building = false;
              }, 1000);
              return { success: true, plugin: plugin.name };
            })
            .catch(error => {
              console.error(`Error building plugin ${plugin.name}:`, error);
              plugin.building = false;
              return { success: false, plugin: plugin.name, error };
            });
        });

        const results = await Promise.all(buildPromises);
        const successful = results.filter(r => r.success).length;
        const failed = results.filter(r => !r.success).length;

        if (failed === 0) {
          alert(`모든 플러그인(${successful}개)이 성공적으로 빌드되었습니다!`);
        } else {
          alert(`${successful}개 플러그인은 성공, ${failed}개 플러그인은 실패했습니다.`);
        }
      } catch (error) {
        console.error('Error building all plugins:', error);
        alert('플러그인 빌드 중 오류가 발생했습니다.');
      }
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
  min-width: 8rem;
  height: 2rem;
  padding: 0.2rem 0.5rem;
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

.setting__button--icon:hover {
  /* 톱니바퀴 이미지 마우스 올렸을 때, 1.1배 커지고 rotate 애니메이션 */
  transition: 0.5s;
  transform: scale(1.1) rotate(90deg);
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
  font-size: 1rem;
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

.build-button {
  padding: 8px 16px;
  margin: 8px 0;
  border: none;
  border-radius: 6px;
  font-size: 0.9rem;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.3s ease;
  background-color: #2196f3;
  color: white;
  min-width: 80px;
}

.build-button:hover:not(:disabled) {
  background-color: #1976d2;
  transform: translateY(-1px);
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.build-button:disabled {
  cursor: not-allowed;
  opacity: 0.6;
}

.build-button.building {
  background-color: #ff9800;
  cursor: not-allowed;
}

.build-button.image-exists {
  background-color: #4caf50;
  cursor: not-allowed;
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

input:checked+.slider.w-color {
  background-color: #ccc;
}

input:checked+.slider.icon {
  background-color: #a37eba;
}

.slider.icon:before {
  background-color: #ffe05d;
}

.slider.icon:after {
  background-color: #e2df23;
}

input:checked+.slider {
  background-color: #2196f3;
}

input:checked+.slider:before {
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

.build-all__button {
  min-width: 8rem;
  height: 2rem;
  padding: 0.2rem 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  background: #2196f3;
  color: white;
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.build-all__button:hover {
  cursor: pointer;
  background: #1976d2;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.build-all__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
  filter: brightness(0) invert(1);
}
</style>
