<template>
  <div v-html="htmlContent"></div>
</template>

<script>
import { getHtml } from "@/api/index";

export default {
  data() {
    return {
      htmlContent: "",
      htmlName: "",
    };
  },
  mounted() {
    this.htmlName = this.$route.query.htmlName;
    this.loadHtmlContent(this.htmlName);
  },
  methods: {
    async loadHtmlContent(htmlName) {
      try {
        const response = await getHtml(htmlName);
        this.htmlContent = response.data;
      } catch (error) {
        console.log(error);
      }
    },
  },
};
</script>
